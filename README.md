# Quantitative Exploration of Elementary Cellular Automata

This repository contains a complete computational study of all 256 Elementary Cellular Automata (ECA) rules using quantitative dynamical metrics.

The project performs a systematic classification based on:

- Mean Shannon entropy
- Normalized sensitivity to initial perturbations
- Lyapunov-like growth rate
- Heuristic detection of shift-periodic localized structures (glider-like behavior)

The goal is to provide a reproducible, metric-based perspective on dynamical regimes in ECA.

## What This Project Does

For each of the 256 Wolfram rules:

- Simulates evolution on a 1D periodic lattice
- Runs multiple independent realizations
- Computes statistical summaries:
  - mean
  - standard deviation
  - 95% confidence interval
- Saves all results to CSV

## Generated Outputs

Running:

```bash
python3 nonlinear.py
````

Produces:

* `eca_metrics.csv` — full statistical dataset
* `entropy_hist.png` — entropy distribution across rules
* `entropy_vs_sensitivity.png` — entropy vs sensitivity scatter plot
* `rules_panel.png` — 2×2 panel of selected space-time diagrams

## Metrics

### Entropy

Average Shannon entropy of configurations over time:

H = -p log2(p) - (1 - p) log2(1 - p)

where p is the fraction of active cells.

### Sensitivity

Normalized mean Hamming distance between two trajectories differing by one bit.

Values lie in [0, 1].

### Lyapunov-like Exponent

λ = (1 / T) * log(d(T) / d(0))

If perturbations fully decay, the value is treated as undefined and excluded from statistics (reported via lyap_bad_frac).

### Glider Heuristic

A simple shift-periodicity test from a single-site initial condition.

This is a coarse indicator and not a formal particle detector.

## Reproducibility

All experiments use fixed random seeds per rule.

Default parameters:

* Lattice size: N = 200
* Time steps: T = 200
* Runs per rule: 20
* Boundary conditions: periodic
* Seed: 12345

These can be modified in nonlinear.py.

## Repository Structure

nonlinear.py                # Main analysis script
eca_metrics.csv             # Output dataset
entropy_hist.png            # Entropy distribution
entropy_vs_sensitivity.png  # Scatter plot
rules_panel.png             # Space-time diagrams
README.md
CITATION.cff

## Requirements

* Python 3.10+
* NumPy
* Matplotlib

Install with:

```bash
pip install numpy matplotlib
```

## Scientific Context

This repository accompanies a quantitative study of ECA dynamical regimes and is intended for reproducible computational experimentation.

It does not attempt to prove universality or physical interpretations; the focus is purely on dynamical classification via statistical metrics.

## License

MIT License.