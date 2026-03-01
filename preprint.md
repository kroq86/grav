# Quantitative Exploration of Elementary Cellular Automata

Kirill Ostapenko  
2026  

---

## Abstract

We perform a systematic computational analysis of all 256 elementary cellular automata (ECA) using quantitative dynamical metrics. For each rule, we compute time-averaged Shannon entropy, normalized perturbation sensitivity, and a Lyapunov-like growth rate. Statistical summaries are obtained from multiple independent realizations with fixed seeds to ensure reproducibility.  

The results provide a fully enumerated metric-based classification of the ECA rule space and reveal quantitative separation between low-entropy, periodic, and high-sensitivity regimes.

---

## 1. Introduction

Elementary cellular automata (ECA), introduced by Stephen Wolfram, consist of one-dimensional binary-state systems with nearest-neighbor interactions and synchronous updates. Despite their minimal definition, they exhibit a wide range of dynamical behaviors including fixed points, periodic structures, chaotic regimes, and complex localized structures.

While qualitative classifications exist, a systematic quantitative survey over the entire 256-rule space provides additional insight into dynamical regimes.

This work performs a full computational enumeration of all ECA rules and evaluates each using reproducible statistical metrics.

---

## 2. Model Definition

An ECA consists of:

- A 1D lattice of size \( N \)
- Binary states \( s_i \in \{0,1\} \)
- Periodic boundary conditions
- Update rule determined by a Wolfram rule number \( r \in [0,255] \)

The update is:

\[
s_i(t+1) = f(s_{i-1}(t), s_i(t), s_{i+1}(t))
\]

where \( f \) is defined by the binary representation of the rule number.

---

## 3. Experimental Setup

For each rule:

- Lattice size: \( N = 200 \)
- Time steps: \( T = 200 \)
- Runs per rule: 20
- Boundary conditions: periodic
- Fixed deterministic seeds

All 256 rules are enumerated exhaustively.

---

## 4. Metrics

### 4.1 Shannon Entropy

For a configuration at time \( t \):

\[
H_t = -p \log_2 p - (1-p)\log_2(1-p)
\]

where \( p \) is the fraction of active cells.

The reported entropy is the time-averaged value across runs.

---

### 4.2 Normalized Sensitivity

Two trajectories differing by a single bit in the initial state are evolved.

Let \( d_t \) be the Hamming distance at time \( t \).

The normalized sensitivity is:

\[
S = \frac{1}{N} \langle d_t \rangle
\]

averaged over time and runs.

---

### 4.3 Lyapunov-like Growth Rate

We define:

\[
\lambda = \frac{1}{T} \log \frac{d(T)}{d(0)}
\]

If perturbations decay completely (\( d(T)=0 \)), the value is excluded from statistical averaging.

---

### 4.4 Glider Heuristic

A shift-periodicity test is applied to single-site initial conditions.  
This heuristic detects translationally propagating structures but is not a formal particle detector.

---

## 5. Results

### 5.1 Entropy Distribution

Entropy values span the full range from near-zero fixed-point rules to high-entropy chaotic rules.

Low-entropy rules correspond to fixed or periodic dynamics.  
High-entropy rules tend to exhibit irregular behavior.

---

### 5.2 Sensitivity and Chaotic Regimes

Rules with high normalized sensitivity display persistent perturbation growth.

Notably:

- Rule 30
- Rule 45
- Rule 86
- Rule 122
- Rule 126

exhibit strong perturbation amplification.

---

### 5.3 Entropy vs Sensitivity Plane

Plotting entropy against sensitivity reveals distinct clusters:

- Low entropy / low sensitivity (ordered rules)
- High entropy / low sensitivity (structured complexity)
- High entropy / high sensitivity (chaotic rules)

This separation provides a quantitative regime classification.

---

### 5.4 Localized Structures

The glider heuristic identifies rules with shift-periodic behavior, including:

- Rule 41
- Rule 97
- Rule 106
- Rule 184

These rules exhibit propagating localized patterns under specific initial conditions.

---

## 6. Discussion

The full enumeration confirms that:

1. Entropy alone does not determine sensitivity.
2. High sensitivity is concentrated in a subset of rules.
3. Localized propagating structures can coexist with moderate entropy.

The metric-based classification aligns partially with classical qualitative classifications but provides quantitative separation.

---

## 7. Reproducibility

All code and datasets are publicly available:

https://github.com/kroq86/grav

Results are deterministic given the fixed seed and parameters.

---

## 8. Conclusion

We provide a complete quantitative survey of all elementary cellular automata using reproducible statistical metrics.

The work establishes:

- A dataset covering all 256 rules
- Statistical characterization of dynamical regimes
- A reproducible computational framework

Future work may extend analysis to larger lattices, longer time horizons, and formal particle detection methods.

---

## References

1. Wolfram, S. "Statistical Mechanics of Cellular Automata." Reviews of Modern Physics (1983).
2. Wolfram, S. *A New Kind of Science* (2002).