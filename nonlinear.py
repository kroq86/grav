from grav.nonlinear_core import *  # noqa: F401,F403


if __name__ == "__main__":
    from grav.nonlinear_core import (
        analyze_all,
        plot_entropy_hist,
        plot_entropy_vs_sensitivity,
        save_csv,
        save_rules_panel,
    )

    results = analyze_all(runs=20, N=200, T=200, seed=12345)
    save_csv(results, "eca_metrics.csv")
    plot_entropy_hist(results, "entropy_hist.png")
    plot_entropy_vs_sensitivity(results, "entropy_vs_sensitivity.png", with_errorbars=True)
    save_rules_panel((30, 41, 106, 184), path="rules_panel.png", N=500, T=500, init="single")
    save_rules_panel((30, 45, 86, 89), path="rules_panel_random.png", N=500, T=500, init="random")
