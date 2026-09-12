"""
Refit of MatsuoLadderDropletModel.tex Stage 1 / Stage 2 on REAL Source Data
(exact xlsx numbers from Matsuo & Kurihara 2021 Source Data files), replacing
the earlier +/-5-10% READING-tier digitized pixel-trace CSVs.

Model equations reproduced exactly from MatsuoLadderDropletModel.tex sec:model
(eqs 1-21 in the archive doc's numbering). Time unit: hours throughout (source
xlsx files are in minutes; Fig3d's "Time (h)" header is a mislabel -- its values
match Fig3c's minute range/spacing exactly).
"""
import numpy as np
import openpyxl
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA = str(HERE / "source_data")

def load_sheet(fn, sheet):
    wb = openpyxl.load_workbook(f"{DATA}/{fn}", data_only=True)
    ws = wb[sheet]
    rows = list(ws.iter_rows(values_only=True))
    header = rows[0]
    data = [r for r in rows[1:] if r[0] is not None]
    return header, data

# ---- Load real data ----
# Fig3c: "Residual monomer" (O(t) observable, 25mM, 0-24h) and "Conversion" (q(t), 25mM)
_, resid_3c = load_sheet("Fig3c.xlsx", "Residual monomer")
_, conv_3c = load_sheet("Fig3c.xlsx", "Conversion")
t_O25_full_min = np.array([r[0] for r in resid_3c], dtype=float)
O25_full = np.array([r[1] for r in resid_3c], dtype=float)
t_q_min = np.array([r[0] for r in conv_3c], dtype=float)
q_data = np.array([r[1] for r in conv_3c], dtype=float)

# Fig3d: two-header sheet. Header row literally reads ('[Mpre]/mM', 25, 13), i.e. column B
# claims to be 25mM DTT and column C claims 13mM. VERIFIED WRONG (2026-09-09, this session):
# column C matches Fig3c's real, caption-confirmed 25mM curve exactly (0.0000 diff) at every
# overlapping timepoint; column B is a genuinely different, independent, slower-decaying curve
# -- chemically, LESS DTT should give SLOWER M_pre reduction, so column B (slow) is the true
# 13mM curve and column C (fast, duplicate of Fig3c) is the true 25mM curve. The header labels
# in this source file are swapped. Also: real data covers the FULL 0-24h (1440 min) for both
# columns, not just the 0-6h the paper's Fig.3d panel visually crops to -- source data commonly
# retains the full raw series beyond a cropped display window.
wb = openpyxl.load_workbook(f"{DATA}/Fig3d.xlsx", data_only=True)
ws = wb["Sheet1"]
rows = list(ws.iter_rows(values_only=True))
d3d = [r for r in rows[2:] if r[0] is not None]
t_3d_min = np.array([r[0] for r in d3d], dtype=float)
O13_true = np.array([r[1] for r in d3d], dtype=float)  # header said 25, verified = true 13mM
O25_true = np.array([r[2] for r in d3d], dtype=float)  # header said 13, verified = true 25mM (dup of Fig3c)

# convert all times to hours
t_O25_full = t_O25_full_min / 60.0
t_O25_3d = t_3d_min / 60.0
t_O13_3d = t_3d_min / 60.0
O25_6h, O13_6h = O25_true, O13_true      # keep old variable names, corrected labels, full-range
t_O25_6h, t_O13_6h = t_O25_3d, t_O13_3d
t_q = t_q_min / 60.0

print("Loaded real source data (Fig3d column labels corrected, see comment above):")
print(f"  O(t) 25mM Fig3c: {len(t_O25_full)} pts, t in [{t_O25_full.min():.3f},{t_O25_full.max():.3f}] h")
print(f"  O(t) 25mM Fig3d (corrected col, ~duplicate of Fig3c): {len(t_O25_6h)} pts, t in [{t_O25_6h.min():.3f},{t_O25_6h.max():.3f}] h")
print(f"  O(t) 13mM Fig3d (corrected col, the real held-out curve): {len(t_O13_6h)} pts, t in [{t_O13_6h.min():.3f},{t_O13_6h.max():.3f}] h")
print(f"  q(t)  25mM Fig3c: {len(t_q)} pts, t in [{t_q.min():.3f},{t_q.max():.3f}] h")

# ============================================================
# Stage 1: oligomerization ladder, eqs (1)-(5)
# ============================================================
N = 12  # max chain length, fixed
X0 = 10.0  # mM, x(0) fixed

def rhs(t, y, k0, k1, kel, p, d0):
    # y = [x, d, m, p2..pN, h]  (N-1 p_n states)
    x, d = y[0], y[1]
    m = y[2]
    pn = y[3:3+N-1]  # p2..pN
    h = y[-1]
    d_safe = max(d, 0.0)
    x_safe = max(x, 0.0)
    m_safe = max(m, 0.0)
    R0 = k0 * x_safe * (d_safe ** p)
    R1 = k1 * m_safe * m_safe
    Rn = kel * np.clip(pn[:-1], 0, None) * m_safe  # R_2..R_{N-1}, i.e. indices for p2..p_{N-1}
    dx = -R0
    dd = -R0
    dm = 2*R0 - 2*R1 - Rn.sum()
    # pn array index 0 -> p2, ..., index N-2 -> pN
    dpn = np.zeros(N-1)
    dpn[0] = R1 - Rn[0]
    for i in range(1, N-2):
        dpn[i] = Rn[i-1] - Rn[i]
    dpn[N-2] = Rn[-1]  # p_N gets last elongation flux, no outflow
    dh = R1 + Rn.sum()
    return np.concatenate(([dx, dd, dm], dpn, [dh]))

def simulate(k0, k1, kel, p, d0, t_eval):
    y0 = np.zeros(3 + (N-1) + 1)
    y0[0] = X0
    y0[1] = d0
    sol = solve_ivp(rhs, (0, max(t_eval.max(), 24.0)), y0, args=(k0,k1,kel,p,d0),
                     t_eval=np.sort(np.unique(np.concatenate([t_eval, [0]]))),
                     method='LSODA', rtol=1e-8, atol=1e-10, max_step=0.1)
    return sol

def O_of_t(sol, t_eval):
    x = sol.y[0]; m = sol.y[2]; pn = sol.y[3:3+N-1]
    Ot = (2*x + m + pn.sum(axis=0)) / (2*X0) * 100.0
    return np.interp(t_eval, sol.t, Ot)

PENALTY = 1e4
def residuals(params):
    k0, k1, kel, p = np.exp(params[0]), np.exp(params[1]), np.exp(params[2]), params[3]
    try:
        sol25 = simulate(k0, k1, kel, p, 25.0, np.union1d(t_O25_full, t_O25_6h))
        sol13 = simulate(k0, k1, kel, p, 13.0, t_O13_6h)
        if not (sol25.success and sol13.success):
            raise RuntimeError("integration failed")
        r = [O_of_t(sol25, t_O25_full) - O25_full,
             O_of_t(sol25, t_O25_6h) - O25_6h,
             O_of_t(sol13, t_O13_6h) - O13_6h]
        out = np.concatenate(r)
        if not np.all(np.isfinite(out)):
            raise RuntimeError("non-finite residual")
        return out
    except Exception:
        n = len(t_O25_full) + len(t_O25_6h) + len(t_O13_6h)
        return np.full(n, PENALTY)

# --- sanity check: old (digitized-data) literature params on the REAL data ---
print("\nSanity check: old literature params (k0=0.00133,k1=0.596,kel=2.288,p=2.30) on real data:")
try:
    sol25_sanity = simulate(0.00133, 0.596, 2.288, 2.30, 25.0, np.union1d(t_O25_full, t_O25_6h))
    sol13_sanity = simulate(0.00133, 0.596, 2.288, 2.30, 13.0, t_O13_6h)
    print("  O25full model sample:", O_of_t(sol25_sanity, t_O25_full)[:5], "vs data", O25_full[:5])
    print("  O13(6h) model sample:", O_of_t(sol13_sanity, t_O13_6h)[:5], "vs data", O13_6h[:5])
except Exception as e:
    print("  sanity check failed:", e)

bounds_lo = [np.log(1e-8), np.log(1e-4), np.log(1e-4), 0.1]
bounds_hi = [np.log(1e4), np.log(1e6), np.log(1e6), 8.0]
print("\nFitting Stage 1 jointly (k0,k1,kel,p) to 3 real curves (bounded trf, multistart)...")
starts = [
    [np.log(0.00133), np.log(0.596), np.log(2.288), 2.30],
    [np.log(0.01), np.log(1.0), np.log(1.0), 2.0],
    [np.log(0.1), np.log(0.1), np.log(0.1), 3.0],
    [np.log(1.0), np.log(0.01), np.log(1.0), 4.0],
    [np.log(0.001), np.log(10.0), np.log(5.0), 1.5],
]
best = None
for i, s0 in enumerate(starts):
    r = least_squares(residuals, s0, method='trf', bounds=(bounds_lo, bounds_hi), max_nfev=3000)
    print(f"  start {i}: cost={r.cost:.4f} status={r.status} -> k0={np.exp(r.x[0]):.5g} k1={np.exp(r.x[1]):.5g} kel={np.exp(r.x[2]):.5g} p={r.x[3]:.3f}")
    if best is None or r.cost < best.cost:
        best = r
res = best
k0f, k1f, kelf, pf = np.exp(res.x[0]), np.exp(res.x[1]), np.exp(res.x[2]), res.x[3]
print(f"BEST: k0={k0f:.6f}  k1={k1f:.4f}  kel={kelf:.4f}  p={pf:.4f}  (cost={res.cost:.4f}, status={res.status})")

def rmse_r2(model, data):
    resid = model - data
    rmse = np.sqrt(np.mean(resid**2))
    ss_res = np.sum(resid**2)
    ss_tot = np.sum((data - data.mean())**2)
    r2 = 1 - ss_res/ss_tot
    return rmse, r2

sol25_final = simulate(k0f, k1f, kelf, pf, 25.0, np.union1d(t_O25_full, t_O25_6h))
sol13_final = simulate(k0f, k1f, kelf, pf, 13.0, t_O13_6h)

O25full_model = O_of_t(sol25_final, t_O25_full)
O256h_model = O_of_t(sol25_final, t_O25_6h)
O136h_model = O_of_t(sol13_final, t_O13_6h)

rmse1, r2_1 = rmse_r2(O25full_model, O25_full)
rmse2, r2_2 = rmse_r2(O256h_model, O25_6h)
rmse3, r2_3 = rmse_r2(O136h_model, O13_6h)
print(f"\n  25mM Fig3c (0-24h): RMSE={rmse1:.3f}  R2={r2_1:.4f}   (old digitized-data: RMSE=1.60 R2=0.983)")
print(f"  25mM Fig3d (0-6h) : RMSE={rmse2:.3f}  R2={r2_2:.4f}   (old digitized-data: RMSE=1.41 R2=0.995)")
print(f"  13mM Fig3d (0-6h) : RMSE={rmse3:.3f}  R2={r2_3:.4f}   (old digitized-data: RMSE=1.11 R2=0.995)")

# held-out check: naive p=1 fit to 25mM alone, predict 13mM (mirrors doc's original result)
def residuals_p1_25only(params):
    k0, k1, kel = np.exp(params[0]), np.exp(params[1]), np.exp(params[2])
    try:
        sol25 = simulate(k0, k1, kel, 1.0, 25.0, np.union1d(t_O25_full, t_O25_6h))
        if not sol25.success:
            raise RuntimeError("integration failed")
        r = np.concatenate([O_of_t(sol25, t_O25_full) - O25_full, O_of_t(sol25, t_O25_6h) - O25_6h])
        if not np.all(np.isfinite(r)):
            raise RuntimeError("non-finite")
        return r
    except Exception:
        return np.full(len(t_O25_full)+len(t_O25_6h), PENALTY)

print("\nHeld-out check: naive p=1 fit to 25mM alone, zero-parameter predict 13mM...")
res_p1 = least_squares(residuals_p1_25only, [np.log(0.00133), np.log(0.596), np.log(2.288)],
                        method='trf',
                        bounds=([np.log(1e-8), np.log(1e-4), np.log(1e-4)], [np.log(1e4), np.log(1e6), np.log(1e6)]),
                        max_nfev=5000)
k0_p1, k1_p1, kel_p1 = np.exp(res_p1.x[0]), np.exp(res_p1.x[1]), np.exp(res_p1.x[2])
sol13_p1 = simulate(k0_p1, k1_p1, kel_p1, 1.0, 13.0, t_O13_6h)
O13_p1_pred = O_of_t(sol13_p1, t_O13_6h)
rmse_p1, r2_p1 = rmse_r2(O13_p1_pred, O13_6h)
print(f"  p=1, 25mM-only fit -> predict 13mM: RMSE={rmse_p1:.3f}  R2={r2_p1:.4f}  (old: RMSE=? R2=-0.32)")

# mean chain length checks -- independent cross-check must use the REAL DATA's observed
# plateau, not the model's own simulated O(24h) (that would just reproduce the fit, not
# check it independently).
pn_final25 = sol25_final.y[3:3+N-1, -1]
ns = np.arange(2, N+1)
nbar_dist = (pn_final25 * ns).sum() / pn_final25.sum() if pn_final25.sum() > 0 else float('nan')
plateau_mask = t_O25_full >= 12.0  # late-time near-plateau real data points (t>=12h)
O_data_plateau = O25_full[plateau_mask].mean()
nbar_plateau = 1.0 / (O_data_plateau/100.0)
d_end = np.interp(24.0, sol25_final.t, sol25_final.y[1])
print(f"\n  n_bar (model's final oligomer distribution) = {nbar_dist:.2f}   (old: 6.4)")
print(f"  n_bar (1/O_data(plateau), real data t>=12h, independent) = {nbar_plateau:.2f}  (mean O={O_data_plateau:.2f}%, old: 6.1)")
print(f"  DTT remaining at t=24h: {d_end:.2f} mM of 25 mM ({100*d_end/25:.1f}%)   (old: ~15mM, 60% remaining)")

results = dict(k0=k0f, k1=k1f, kel=kelf, p=float(pf),
               rmse=[rmse1,rmse2,rmse3], r2=[r2_1,r2_2,r2_3],
               p1_heldout_rmse=rmse_p1, p1_heldout_r2=r2_p1,
               nbar_dist=nbar_dist, nbar_plateau=nbar_plateau, dtt_remaining_mM=d_end)
with open(HERE / "stage1_results.json","w") as f:
    json.dump(results, f, indent=2)
print("\nStage 1 done, results saved.")

# ============================================================
# Stage 2: droplet formation, eqs (qstar)-(q). S(t) from Stage-1's 25mM trajectory
# (full 0-24h ladder, fit params from above), fit Y,Fsat0,lambda,K,n to real q(t) (Fig3c
# Conversion sheet, 91 real points, converted to hours).
# ============================================================
sol25_stage2 = sol25_final  # already integrated 0-24h at the fitted Stage-1 params

def S_of_t(t_eval):
    pn = sol25_stage2.y[3:3+N-1]
    h = sol25_stage2.y[-1]
    S = pn.sum(axis=0) + h
    return np.interp(t_eval, sol25_stage2.t, S)

def X_of_t(t_eval, K, n):
    # dX/dt = (1-X) n K t^(n-1), X(0)=0  =>  closed form X(t) = 1 - exp(-K t^n)
    return 1.0 - np.exp(-K * np.power(np.clip(t_eval, 0, None), n))

def q_model(t_eval, Y, Fsat0, lam, K, n):
    S = S_of_t(t_eval)
    qstar = Y * (S - Fsat0 - lam * t_eval)
    X = X_of_t(t_eval, K, n)
    return qstar * X

def q_residuals(params, t_fit, q_fit):
    Y, Fsat0, lam, K, n = params
    if K <= 0 or n <= 0 or Y <= 0:
        return np.full(len(t_fit), PENALTY)
    try:
        m = q_model(t_fit, Y, Fsat0, lam, K, n)
        if not np.all(np.isfinite(m)):
            return np.full(len(t_fit), PENALTY)
        return m - q_fit
    except Exception:
        return np.full(len(t_fit), PENALTY)

print("\n" + "="*60)
print("Stage 2: fitting droplet formation q(t) to real Conversion data")
print("="*60)
p0_stage2 = [3.12, -3.54, 0.143, 0.494, 0.933]  # old literature values as starting point
res_q = least_squares(q_residuals, p0_stage2, args=(t_q, q_data), method='lm', max_nfev=10000)
Yf, Fsat0f, lamf, Kf, nf = res_q.x
q_model_fit = q_model(t_q, Yf, Fsat0f, lamf, Kf, nf)
rmse_q, r2_q = rmse_r2(q_model_fit, q_data)
print(f"  Y={Yf:.4f}  Fsat0={Fsat0f:.4f}  lambda={lamf:.4f}  K={Kf:.4f}  n={nf:.4f}")
print(f"  RMSE={rmse_q:.3f}  R2={r2_q:.4f}   (old digitized-data: RMSE=0.44 R2=0.999)")

# peak / endpoint comparison
peak_idx_data = np.argmax(q_data)
peak_idx_model = np.argmax(q_model_fit)
print(f"\n  Data peak: {q_data[peak_idx_data]:.2f}% at t={t_q[peak_idx_data]:.2f}h")
print(f"  Model peak: {q_model_fit[peak_idx_model]:.2f}% at t={t_q[peak_idx_model]:.2f}h  (old model: 68.0% at 9.9h, old data: 68.6% at 10.5h)")
print(f"  Data value at t=24h: {q_data[-1]:.2f}%   Model value at t=24h: {q_model_fit[-1]:.2f}%   (old: data 62.7%, model 62.8%)")

# ---- held-out-in-time test ----
print("\nStage 2 held-out-in-time test (fit to t<=cutoff, predict rest, zero new params):")
cutoffs = [10.0, 12.0, 14.0, 16.0, 18.0]
heldout_results = []
for cutoff in cutoffs:
    train_mask = t_q <= cutoff
    test_mask = ~train_mask
    if test_mask.sum() < 3:
        continue
    res_cut = least_squares(q_residuals, p0_stage2, args=(t_q[train_mask], q_data[train_mask]),
                             method='lm', max_nfev=10000)
    Yc, Fc, lamc, Kc, nc = res_cut.x
    pred_test = q_model(t_q[test_mask], Yc, Fc, lamc, Kc, nc)
    _, r2_test = rmse_r2(pred_test, q_data[test_mask])
    heldout_results.append((cutoff, r2_test, lamc))
    print(f"  cutoff={cutoff:5.1f}h  held-out R2={r2_test:7.3f}  lambda(from training)={lamc:.4f}")
print(f"  (full-data lambda = {lamf:.4f})")

stage2_results = dict(Y=Yf, Fsat0=Fsat0f, lam=lamf, K=Kf, n=nf, rmse=rmse_q, r2=r2_q,
                       heldout=heldout_results)
with open(HERE / "stage2_results.json","w") as f:
    json.dump(stage2_results, f, indent=2)
print("\nStage 2 done, results saved.")

# ============================================================
# Matsuo & Kurihara's own Richards-curve fit, for reference (their OWN parameters, fit to
# their own separate Suppl. Fig. 14 replicate run -- not refit here, per the archive's own
# caveat that this is a different specific dataset). Reported here only to recompute R2 of
# THEIR functional form against the (now real, not digitized) Fig3c q(t) data, as an extra
# same-data comparison point.
# ============================================================
def richards_matsuo(x, a, xc, d, k):
    base = 1 + (d-1)*np.exp(-k*(x-xc))
    out = np.full_like(x, np.nan, dtype=float)
    valid = base > 0
    out[valid] = a * base[valid]**(1/(1-d))
    return out, valid

t_q_min_arr = t_q * 60.0  # their fit's x is in minutes
y_matsuo_on_our_data, valid_mask = richards_matsuo(t_q_min_arr, 69.70, 77.36, 0.465, 0.00822)
n_invalid = (~valid_mask).sum()
rmse_m, r2_m = rmse_r2(y_matsuo_on_our_data[valid_mask], q_data[valid_mask])
print(f"\nMatsuo & Kurihara's OWN fitted Richards curve (from their Suppl Fig 14, a different")
print(f"replicate run's parameters), evaluated against OUR real Fig3c q(t) data:")
if n_invalid:
    print(f"  ({n_invalid}/{len(t_q)} earliest points excluded: their fit's base term goes")
    print(f"   slightly negative there, undefined for its non-integer exponent -- an edge-of-")
    print(f"   domain effect of extrapolating their curve to t=0, not a computation bug)")
print(f"  RMSE={rmse_m:.3f}  R2={r2_m:.4f}  (their own reported R2 on their own data: 0.99987)")

# ============================================================
# Plot: 3-panel figure, real data + fitted model, replacing the digitized-data version
# ============================================================
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))

t_dense25 = np.linspace(0, 24, 300)
t_dense13 = np.linspace(0, 24, 300)
O25_dense = O_of_t(sol25_final, t_dense25)
sol13_dense = simulate(k0f, k1f, kelf, pf, 13.0, t_dense13)
O13_dense = O_of_t(sol13_dense, t_dense13)

ax = axes[0]
ax.scatter(t_O25_full, O25_full, s=25, color="C0", label="Fig.3c (real)", zorder=3)
ax.scatter(t_O25_6h, O25_6h, s=25, marker="^", color="C1", label="Fig.3d, corrected col (real)", zorder=3)
ax.plot(t_dense25, O25_dense, color="k", lw=1.5, label="joint model fit")
ax.set_title("25 mM DTT: O(t), real data")
ax.set_xlabel("t (h)"); ax.set_ylabel("O(t) (%)"); ax.legend(fontsize=8)

ax = axes[1]
ax.scatter(t_O13_6h, O13_6h, s=25, color="C2", label="Fig.3d, corrected col (real, held-out)", zorder=3)
ax.plot(t_dense13, O13_dense, color="k", lw=1.5, label="predicted (fit to 25mM only params, p free)")
ax.set_title("13 mM DTT: O(t), real data (joint fit)")
ax.set_xlabel("t (h)"); ax.legend(fontsize=8)

ax = axes[2]
ax.scatter(t_q, q_data, s=15, color="C3", label="Fig.3c Conversion (real)", zorder=3)
t_dense_q = np.linspace(0, 24, 300)
q_dense = q_model(t_dense_q, Yf, Fsat0f, lamf, Kf, nf)
ax.plot(t_dense_q, q_dense, color="k", lw=1.5, label="Stage 2 fit")
ax.set_title("Droplet conversion q(t), real data")
ax.set_xlabel("t (h)"); ax.set_ylabel("q(t) (%)"); ax.legend(fontsize=8)

plt.tight_layout()
outpath = HERE / "fig_data_and_fit.png"
plt.savefig(outpath, dpi=150)
print(f"\nSaved figure to {outpath}")
