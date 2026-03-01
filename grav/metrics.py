from __future__ import annotations

import argparse
from pathlib import Path

from grav.nonlinear_core import (
    analyze_all,
    plot_entropy_hist,
    plot_entropy_vs_sensitivity,
    save_csv,
    save_rules_panel,
)


def _parse_rules(value: str) -> tuple[int, int, int, int]:
    parts = tuple(int(item.strip()) for item in value.split(","))
    if len(parts) != 4:
        raise argparse.ArgumentTypeError("expected exactly 4 comma-separated rule numbers")
    return parts


def main() -> None:
    parser = argparse.ArgumentParser(description="Compute ECA metrics and save reproducible outputs.")
    parser.add_argument("--runs", type=int, default=20)
    parser.add_argument("--N", type=int, default=200)
    parser.add_argument("--T", type=int, default=200)
    parser.add_argument("--seed", type=int, default=12345)
    parser.add_argument("--out-dir", default="results")
    parser.add_argument("--csv", default="eca_metrics.csv")
    parser.add_argument("--entropy-hist", default="entropy_hist.png")
    parser.add_argument("--entropy-vs-sensitivity", default="entropy_vs_sensitivity.png")
    parser.add_argument("--rules-panel", default="rules_panel.png")
    parser.add_argument("--rules-panel-random", default="rules_panel_random.png")
    parser.add_argument("--rules", type=_parse_rules, default=(30, 41, 106, 184))
    parser.add_argument("--random-rules", type=_parse_rules, default=(30, 45, 86, 89))
    parser.add_argument("--skip-plots", action="store_true")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    results = analyze_all(runs=args.runs, N=args.N, T=args.T, seed=args.seed)
    save_csv(results, str(out_dir / args.csv))

    if args.skip_plots:
        return

    plot_entropy_hist(results, str(out_dir / args.entropy_hist))
    plot_entropy_vs_sensitivity(results, str(out_dir / args.entropy_vs_sensitivity), with_errorbars=True)
    save_rules_panel(args.rules, path=str(out_dir / args.rules_panel), N=500, T=500, init="single")
    save_rules_panel(
        args.random_rules,
        path=str(out_dir / args.rules_panel_random),
        N=500,
        T=500,
        init="random",
    )


if __name__ == "__main__":
    main()
