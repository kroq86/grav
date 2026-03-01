from __future__ import annotations

import argparse
from pathlib import Path

import grav.cluster_core as legacy_cluster


def main() -> None:
    parser = argparse.ArgumentParser(description="Cluster ECA metrics from a CSV file.")
    parser.add_argument("--csv", default="results/eca_metrics.csv")
    parser.add_argument("--out-dir", default="results")
    parser.add_argument("--k-min", type=int, default=2)
    parser.add_argument("--k-max", type=int, default=8)
    parser.add_argument("--seed", type=int, default=12345)
    parser.add_argument("--perm-B", type=int, default=500)
    parser.add_argument("--high-entropy-only", action="store_true")
    parser.add_argument("--high-entropy-threshold", type=float, default=0.80)
    parser.add_argument("--features", nargs="+", default=["entropy_mean", "sens_mean"])
    args = parser.parse_args()

    cfg = legacy_cluster.Config(
        csv_path=args.csv,
        out_prefix="",
        features=tuple(args.features),
        k_min=args.k_min,
        k_max=args.k_max,
        seed=args.seed,
        stability_seeds=(0, 1, 2, 3, 4, 5, 6, 7, 8, 9),
        perm_B=args.perm_B,
        high_entropy_only=bool(args.high_entropy_only),
        high_entropy_threshold=float(args.high_entropy_threshold),
    )

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    df0 = legacy_cluster._load_df(cfg.csv_path)
    df = legacy_cluster._select(df0, cfg)

    if cfg.high_entropy_only:
        print(f"High-entropy rules: {len(df)}")

    X = legacy_cluster._make_X(df, cfg)
    scores = legacy_cluster._silhouette_by_k(X, seed=cfg.seed, k_min=cfg.k_min, k_max=cfg.k_max)
    print(f"Silhouette scores: {scores}")

    k_star = legacy_cluster._best_k(scores)
    print(f"Selected k* = {k_star}")

    labels = legacy_cluster._fit_kmeans(X, k_star, cfg.seed)
    ari_mean, ari_std, ari_list = legacy_cluster._stability_ari(X, k_star, cfg.stability_seeds)
    print(f"Stability over seeds (ARI): mean={ari_mean:.3f}, std={ari_std:.3f}")
    if ari_list:
        print(f"ARI list: {ari_list}")

    real_sil = float(legacy_cluster.silhouette_score(X, labels))
    print(f"Silhouette real={real_sil:.3f}")

    perm, p = legacy_cluster._perm_test_silhouette(X, k_star, real_sil, cfg.perm_B, cfg.seed + 999)
    print(f"Permutation p-value={p:.4f}")

    prefix = "refined_" if cfg.high_entropy_only else ""
    out_csv = out_dir / f"{prefix}clusters.csv"
    out_scatter = out_dir / f"{prefix}clusters_HS.png"
    out_perm = out_dir / f"{prefix}permtest.png"

    df_out = df.copy()
    df_out["cluster"] = labels.astype(int)
    df_out.to_csv(out_csv, index=False)

    legacy_cluster._plot_clusters(df_out, X, labels, cfg, str(out_scatter))
    legacy_cluster._plot_perm_hist(perm, real_sil, p, str(out_perm))

    print(f"Saved: {out_scatter}")
    print(f"Saved: {out_csv}")
    print(f"Saved: {out_perm}")


if __name__ == "__main__":
    main()
