# Research Direction

## Scope

This repository is a computational and methodological base for a broader research program on elementary cellular automata (ECA), with two connected goals:

1. reproducible quantitative exploration of the full 256-rule ECA space;
2. formal study of localized structures, domains, interfaces, and translating defects in selected rules.

The repository should be read as an evolving research artifact, not as a completed theory.

## What The Repository Already Supports

### 1. Reproducible metric pipeline for all 256 ECA rules

The current codebase provides:

- entropy-based observables,
- sensitivity to one-bit perturbations,
- a Lyapunov-like growth heuristic,
- a coarse glider-like indicator,
- clustering in feature space,
- stability analysis via adjusted Rand index,
- permutation testing for silhouette significance.

These components are implemented as a reproducible pipeline through:

- `python3 -m grav.metrics --out-dir results`
- `python3 -m grav.cluster --csv results/eca_metrics.csv --out-dir results`

This part of the repository supports a quantitative, data-driven map of ECA rule space. It does not by itself constitute a formal classification theorem.

### 2. Verified computational core

The repository includes a correctness check aligning:

- a list-based ring implementation of ECA evolution,
- a bitwise implementation used in the faster search code.

This ensures that exploratory searches for fronts and defects are grounded in a consistent evolution rule.

### 3. Exploratory search for localized structures in selected rules

The research scripts under `grav/research/` support bounded search for:

- periodic domains of small period,
- finite defects over two-domain backgrounds,
- exact translating fronts satisfying one-step shift conditions,
- candidate particle-like and interface-like structures in selected rules.

At present, this layer is exploratory and rule-specific. It is strongest for the Rule 69 line of investigation.

## What The Repository Does Not Yet Prove

The current repository does not yet provide:

- a complete formal theory of domains in ECA,
- a general particle classification theorem across ECA,
- a complete collision algebra for localized structures,
- a proof that metric clusters correspond to structural classes,
- a finished symbolic-dynamics treatment of interfaces and defects.

These are research goals, not completed results.

## Research Thesis

The long-term hypothesis is that ECA can serve as a minimal laboratory for rigorous study of emergence in discrete nonlinear systems. In particular:

- domains can be defined as spatially or space-time periodic backgrounds,
- particles can be defined as finite defects over such backgrounds satisfying shift-periodicity conditions of the form `F^p(x) = sigma^v(x)`,
- rule-space metrics may correlate with the presence of domains, interfaces, and stable localized structures.

The repository currently supports the computational side of this thesis and partially supports the structural side through bounded search tools.

## Planned Research Trajectory

### Stage 1. Computationally reproducible exploration

- finalize stable metric computation across all 256 rules;
- maintain reproducible outputs and minimal dependencies;
- keep implementation tests aligned with the evolution model.

### Stage 2. Formalization of domains and particles

- define domains as periodic or shift-periodic backgrounds;
- define finite defects over domains;
- formalize particles and oscillators via exact shift-periodicity;
- replace heuristic glider language with precise dynamical definitions.

### Stage 3. Structural analysis in selected rules

- classify periodic domains for selected nontrivial rules;
- search for admissible interface defects;
- measure periods, velocities, and locality constraints;
- analyze candidate collisions and persistence mechanisms.

### Stage 4. Link structure and metrics

- compare rigorously detected structures with metric observables;
- test whether structural richness aligns with regions of metric rule space;
- determine which observables are genuinely informative and which are only heuristic.

## Positioning

This project is not framed as new fundamental physics. Its scientifically defensible scope is:

- symbolic dynamics,
- discrete nonlinear systems,
- complex systems,
- computational exploration of emergence,
- quantitative classification of cellular automata.

If developed rigorously, the natural research outcome is not a claim about gravity or cosmology, but a theory of localized structures and regime classification in discrete dynamical media, with ECA as the minimal testbed.
