# A Mechanistic Kinetic Model of Peptide Oligomerization and Droplet Formation in the Matsuo–Kurihara System

Alex Shoshitaishvili, Independent Researcher

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22717387.svg)](https://doi.org/10.5281/zenodo.22717387)

## Summary

Matsuo and Kurihara's peptide-coacervate system converts a disulfide precursor
M<sub>pre</sub> into monomer M, which oligomerizes into short peptide chains, which
phase-separate into liquid droplets. This note gives a model that connects those two
stages: a mechanistic reaction-network model of the oligomerization step, feeding a
phase-separation model of droplet formation whose response climbs toward completion
gradually rather than instantaneously.

All fits use the paper's own Source Data files. Fit jointly to all three available real
curves — two DTT concentrations, two figures — the oligomerization stage reproduces every
one of them with a single shared set of rate constants (R²=0.990–0.996); a separate,
genuinely held-out test (fitting the 25 mM data alone under the naive assumption that DTT
enters at first order, then predicting the 13 mM curve with zero new parameters) gives
R²=0.61, which the joint fit resolves by identifying a DTT reaction order of p≈2.37
instead. The droplet-formation stage reproduces the full rise-peak-decline shape of the
droplet-conversion curve essentially exactly (R²=0.9997) — but the paper is explicit that
this second stage has not been shown to have real predictive power with the data
currently available, and states precisely what data would resolve that.

## Files

- `MatsuoLadderDropletModel.tex` — LaTeX source.
- `MatsuoLadderDropletModel.pdf` — compiled paper.
- `fig_data_and_fit.png` — the paper's data/fit figure.
- `matsuo_refit.py` — the fitting script, for reproducibility (loads the source data,
  fits both stages, runs the held-out tests, regenerates the figure).

## Building

```
pdflatex MatsuoLadderDropletModel.tex
pdflatex MatsuoLadderDropletModel.tex
```

(A second pass resolves cross-references and citations.)

## Citing this work

See [`CITATION.cff`](CITATION.cff), or cite via the Zenodo concept DOI (always resolves to the
latest version): [10.5281/zenodo.22717387](https://doi.org/10.5281/zenodo.22717387).

## License

Released under [CC BY 4.0](LICENSE) — you may share and adapt this work for any purpose,
provided you give appropriate credit.
