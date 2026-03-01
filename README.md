Quantitative Exploration of Elementary Cellular Automata

This repository contains a complete computational study of all 256 Elementary Cellular Automata (ECA) rules using quantitative dynamical metrics.

The project performs a systematic classification based on:

Mean Shannon entropy

Normalized sensitivity to initial perturbations

Lyapunov-like growth rate

Heuristic detection of shift-periodic localized structures (glider-like behavior)

The goal is to provide a reproducible, metric-based perspective on dynamical regimes in ECA.

What This Project Does

For each of the 256 Wolfram rules:

Simulates evolution on a 1D periodic lattice.

Runs multiple independent realizations.

Computes statistical summaries:

mean

standard deviation

95% confidence interval

Saves all results to CSV.

Generated Outputs

Running:

python3 nonlinear.py

Produces:

eca_metrics.csv — full statistical dataset

entropy_hist.png — entropy distribution across rules

entropy_vs_sensitivity.png — entropy vs sensitivity scatter plot

rules_panel.png — 2×2 panel of selected space-time diagrams

Metrics
Entropy

Average Shannon entropy of configurations over time:

Sensitivity

Normalized mean Hamming distance between two trajectories differing by one bit.

Values lie in [0,1].

Lyapunov-like Exponent
	​


If perturbations fully decay, the value is treated as undefined and excluded from statistics (reported via lyap_bad_frac).

Glider Heuristic

A simple shift-periodicity test from a single-site initial condition.
This is a coarse indicator and not a formal particle detector.

Reproducibility

All experiments use fixed random seeds per rule.
Default parameters:

Lattice size: N = 200

Time steps: T = 200

Runs per rule: 20

These can be modified in nonlinear.py.

Repository Structure
nonlinear.py        # Main analysis script
eca_metrics.csv     # Output dataset
entropy_hist.png    # Entropy distribution
entropy_vs_sensitivity.png
rules_panel.png
Requirements

Python 3.10+

NumPy

Matplotlib

Install with:

pip install numpy matplotlib

Scientific Context

This repository accompanies a quantitative study of ECA dynamical regimes and is intended for reproducible computational experimentation.

It does not attempt to prove universality or physical interpretations; the focus is purely on dynamical classification via statistical metrics.

License

MIT License.

## Experimental Setup

- Lattice size: N = 200
- Time steps: T = 200
- Runs per rule: 20
- Boundary conditions: periodic
- Seed: 12345