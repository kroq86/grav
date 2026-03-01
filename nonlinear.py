import csv
import numpy as np
import matplotlib.pyplot as plt

# ============================================================
# 1) ECA: rule + simulation
# ============================================================

def make_rule(rule_number: int):
    """
    Wolfram ECA rule (0..255), radius=1, binary states.
    Periodic boundary conditions via np.roll.
    """
    table = np.array([int(x) for x in np.binary_repr(rule_number, width=8)], dtype=np.int8)
    table = table[::-1]  # Wolfram convention

    def step(state: np.ndarray) -> np.ndarray:
        left = np.roll(state, 1)
        right = np.roll(state, -1)
        idx = (left << 2) | (state << 1) | right
        return table[idx]

    return step


def run(step_func, initial_state: np.ndarray, steps: int) -> np.ndarray:
    """
    Run CA for `steps` updates, returns history shape (steps+1, N).
    """
    history = np.empty((steps + 1, initial_state.size), dtype=np.int8)
    state = initial_state.astype(np.int8, copy=True)
    history[0] = state
    for t in range(steps):
        state = step_func(state)
        history[t + 1] = state
    return history


# ============================================================
# 2) Metrics
# ============================================================

def entropy_of_state(state: np.ndarray) -> float:
    p = float(np.mean(state))
    if p <= 0.0 or p >= 1.0:
        return 0.0
    return float(-(p * np.log2(p) + (1.0 - p) * np.log2(1.0 - p)))


def average_entropy(history: np.ndarray) -> float:
    return float(np.mean([entropy_of_state(row) for row in history]))


def sensitivity_normalized(step_func, rng: np.random.Generator, N: int = 200, T: int = 200) -> float:
    """
    Mean Hamming distance over time between two trajectories
    starting from states that differ by 1 bit, normalized by N.
    Returns value in [0, 1].
    """
    s1 = rng.integers(0, 2, size=N, dtype=np.int8)
    s2 = s1.copy()
    s2[rng.integers(0, N)] ^= 1

    h1 = run(step_func, s1, T)
    h2 = run(step_func, s2, T)

    distances = np.sum(h1 != h2, axis=1)  # length T+1
    return float(np.mean(distances) / N)


def lyapunov_like(step_func, rng: np.random.Generator, N: int = 200, T: int = 200) -> float:
    """
    A Lyapunov-like growth rate estimate:
      lambda = (1/T) * log( d(T) / d(0) ), with d(0)=1.
    If d(T)=0 (perturbation dies out), returns -inf.
    """
    s1 = rng.integers(0, 2, size=N, dtype=np.int8)
    s2 = s1.copy()
    s2[rng.integers(0, N)] ^= 1

    h1 = run(step_func, s1, T)
    h2 = run(step_func, s2, T)

    d0 = 1.0
    dT = float(np.sum(h1[-1] != h2[-1]))
    if dT == 0.0:
        return float("-inf")
    return float((1.0 / T) * np.log(dT / d0))


def find_glider_heuristic(step_func, width: int = 41, steps: int = 60, max_shift: int = 10) -> bool:
    """
    Very coarse heuristic:
    start from single 1 in the center; if at some time the entire row equals
    a shifted version of the initial row, flag as glider-like.
    (This is NOT a rigorous glider detector; it's a quick indicator.)
    """
    center = width // 2
    state = np.zeros(width, dtype=np.int8)
    state[center] = 1

    history = run(step_func, state, steps)

    for t in range(2, steps + 1):
        for shift in range(-max_shift, max_shift + 1):
            if shift == 0:
                continue
            if np.array_equal(history[t], np.roll(history[0], shift)):
                return True
    return False


# ============================================================
# 3) Robust statistics
# ============================================================

def mean_std_ci95(x: np.ndarray):
    """
    Returns (mean, std, ci95, bad_frac), where bad are non-finite values.
    Useful for lyapunov_like which can be -inf.
    """
    x = np.asarray(x, dtype=float)
    finite = np.isfinite(x)
    xf = x[finite]
    bad_frac = 1.0 - (xf.size / x.size)

    if xf.size == 0:
        return float("nan"), float("nan"), float("nan"), bad_frac

    mean = float(np.mean(xf))
    std = float(np.std(xf, ddof=1)) if xf.size > 1 else 0.0
    ci = float(1.96 * std / np.sqrt(xf.size)) if xf.size > 1 else 0.0
    return mean, std, ci, bad_frac


# ============================================================
# 4) Analysis
# ============================================================

def analyze_rule(rule: int, runs: int = 20, N: int = 200, T: int = 200, seed: int = 12345):
    step_func = make_rule(rule)

    # Reproducible but different per rule
    rng = np.random.default_rng(seed + rule * 1000003)

    ent = np.empty(runs, dtype=float)
    sens = np.empty(runs, dtype=float)
    lyap = np.empty(runs, dtype=float)

    for k in range(runs):
        init = rng.integers(0, 2, size=N, dtype=np.int8)
        hist = run(step_func, init, T)

        ent[k] = average_entropy(hist)
        sens[k] = sensitivity_normalized(step_func, rng, N=N, T=T)
        lyap[k] = lyapunov_like(step_func, rng, N=N, T=T)

    glider = find_glider_heuristic(step_func)

    ent_m, ent_sd, ent_ci, ent_bad = mean_std_ci95(ent)
    s_m, s_sd, s_ci, s_bad = mean_std_ci95(sens)
    l_m, l_sd, l_ci, l_bad = mean_std_ci95(lyap)

    return {
        "rule": rule,
        "entropy_mean": ent_m, "entropy_std": ent_sd, "entropy_ci95": ent_ci,
        "sens_mean": s_m, "sens_std": s_sd, "sens_ci95": s_ci,
        "lyap_mean": l_m, "lyap_std": l_sd, "lyap_ci95": l_ci, "lyap_bad_frac": l_bad,
        "glider_heuristic": int(glider),
        "runs": runs, "N": N, "T": T, "seed": seed
    }


def analyze_all(runs: int = 20, N: int = 200, T: int = 200, seed: int = 12345):
    results = []
    for rule in range(256):
        r = analyze_rule(rule, runs=runs, N=N, T=T, seed=seed)
        results.append(r)

        print(
            f"Rule {rule:3d}: "
            f"H={r['entropy_mean']:.3f}±{r['entropy_ci95']:.3f}, "
            f"S={r['sens_mean']:.3f}±{r['sens_ci95']:.3f}, "
            f"L={r['lyap_mean']:.3g}±{r['lyap_ci95']:.2g} (bad={r['lyap_bad_frac']:.0%}), "
            f"glider={bool(r['glider_heuristic'])}"
        )

    return results


# ============================================================
# 5) Output: CSV + plots
# ============================================================

def save_csv(results, path="eca_metrics.csv"):
    fieldnames = list(results[0].keys())
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(results)
    print(f"Saved: {path}")


def plot_entropy_hist(results, path="entropy_hist.png"):
    H = np.array([r["entropy_mean"] for r in results], dtype=float)

    plt.figure()
    plt.hist(H, bins=30)
    plt.xlabel("Mean entropy")
    plt.ylabel("Number of rules")
    plt.title("ECA: distribution of mean entropy (over runs)")
    plt.tight_layout()
    plt.savefig(path, dpi=200)
    plt.close()
    print(f"Saved: {path}")


def plot_entropy_vs_sensitivity(results, path="entropy_vs_sensitivity.png", with_errorbars=True):
    H = np.array([r["entropy_mean"] for r in results], dtype=float)
    S = np.array([r["sens_mean"] for r in results], dtype=float)
    Hci = np.array([r["entropy_ci95"] for r in results], dtype=float)
    Sci = np.array([r["sens_ci95"] for r in results], dtype=float)
    G = np.array([r["glider_heuristic"] for r in results], dtype=int)

    plt.figure()
    idx0 = (G == 0)
    idx1 = (G == 1)

    if with_errorbars:
        plt.errorbar(H[idx0], S[idx0], xerr=Hci[idx0], yerr=Sci[idx0],
                     fmt="o", markersize=3, capsize=2, linestyle="none",
                     label="no glider (heur.)")
        plt.errorbar(H[idx1], S[idx1], xerr=Hci[idx1], yerr=Sci[idx1],
                     fmt="^", markersize=4, capsize=2, linestyle="none",
                     label="glider (heur.)")
    else:
        plt.scatter(H[idx0], S[idx0], s=12, marker="o", label="no glider (heur.)")
        plt.scatter(H[idx1], S[idx1], s=18, marker="^", label="glider (heur.)")

    plt.xlabel("Mean entropy")
    plt.ylabel("Mean normalized sensitivity")
    plt.title("ECA: entropy vs sensitivity (means ± CI95)")
    plt.legend()
    plt.tight_layout()
    plt.savefig(path, dpi=200)
    plt.close()
    print(f"Saved: {path}")

def render_rule_history(rule: int, N: int = 400, T: int = 400, init: str = "single", seed: int = 12345) -> np.ndarray:
    """
    Returns history (T+1, N) for a given rule using a chosen initial condition.
    init:
      - "single": single 1 in the center
      - "random": random initial state
    """
    step_func = make_rule(rule)

    if init == "single":
        state0 = np.zeros(N, dtype=np.int8)
        state0[N // 2] = 1
    elif init == "random":
        rng = np.random.default_rng(seed + rule * 1009)
        state0 = rng.integers(0, 2, size=N, dtype=np.int8)
    else:
        raise ValueError("init must be 'single' or 'random'")

    return run(step_func, state0, T)


def save_rule_image(rule: int, path: str, N: int = 400, T: int = 400, init: str = "single"):
    hist = render_rule_history(rule, N=N, T=T, init=init)

    plt.figure(figsize=(6, 6))
    plt.imshow(hist, cmap="binary", interpolation="nearest", aspect="auto")
    plt.title(f"Rule {rule} ({init})")
    plt.xlabel("Cell index")
    plt.ylabel("Time step")
    plt.tight_layout()
    plt.savefig(path, dpi=250)
    plt.close()
    print(f"Saved: {path}")


def save_rules_panel(rules=(30, 41, 106, 184), path: str = "rules_panel.png",
                     N: int = 500, T: int = 500, init: str = "single"):
    """
    Save a 2x2 panel image for a given tuple/list of 4 rules.
    """
    if len(rules) != 4:
        raise ValueError("rules must contain exactly 4 rule numbers for a 2x2 panel")

    fig, axes = plt.subplots(2, 2, figsize=(10, 10))
    axes = axes.ravel()

    for ax, rule in zip(axes, rules):
        hist = render_rule_history(rule, N=N, T=T, init=init)
        ax.imshow(hist, cmap="binary", interpolation="nearest", aspect="auto")
        ax.set_title(f"Rule {rule}")
        ax.set_xlabel("i")
        ax.set_ylabel("t")

    fig.suptitle(f"ECA space-time diagrams ({init} initial condition)", y=0.98)
    fig.tight_layout()
    fig.savefig(path, dpi=250)
    plt.close(fig)
    print(f"Saved: {path}")

# ============================================================
# 6) Main
# ============================================================

if __name__ == "__main__":
    results = analyze_all(runs=20, N=200, T=200, seed=12345)
    save_csv(results, "eca_metrics.csv")
    plot_entropy_hist(results, "entropy_hist.png")
    plot_entropy_vs_sensitivity(results, "entropy_vs_sensitivity.png", with_errorbars=True)
    save_rules_panel((30, 41, 106, 184), path="rules_panel.png", N=500, T=500, init="single")
    save_rules_panel((30, 45, 86, 89), path="rules_panel_random.png", N=500, T=500, init="random")