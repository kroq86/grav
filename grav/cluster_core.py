#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import math
from dataclasses import dataclass

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from sklearn.metrics.cluster import adjusted_rand_score
from sklearn.preprocessing import StandardScaler


@dataclass(frozen=True)
class Config:
    csv_path: str
    out_prefix: str
    features: tuple[str, ...]
    k_min: int
    k_max: int
    seed: int
    stability_seeds: tuple[int, ...]
    perm_B: int
    high_entropy_only: bool
    high_entropy_threshold: float


def _load_df(path: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    if "rule" not in df.columns:
        raise ValueError("CSV должен содержать колонку 'rule'.")
    return df


def _select(df: pd.DataFrame, cfg: Config) -> pd.DataFrame:
    missing = [c for c in cfg.features if c not in df.columns]
    if missing:
        raise ValueError(f"В CSV нет колонок: {missing}")

    out = df.copy()

    if cfg.high_entropy_only:
        if "entropy_mean" not in out.columns:
            raise ValueError("Для high-entropy фильтра нужна колонка 'entropy_mean'.")
        out = out[out["entropy_mean"] >= cfg.high_entropy_threshold].copy()

    out = out.dropna(subset=list(cfg.features)).copy()
    out = out.sort_values("rule").reset_index(drop=True)
    return out


def _make_X(df: pd.DataFrame, cfg: Config) -> np.ndarray:
    X = df.loc[:, list(cfg.features)].to_numpy(dtype=float)
    X = StandardScaler().fit_transform(X)
    return X


def _silhouette_by_k(X: np.ndarray, seed: int, k_min: int, k_max: int) -> dict[int, float]:
    scores: dict[int, float] = {}
    for k in range(k_min, k_max + 1):
        km = KMeans(n_clusters=k, n_init=50, random_state=seed)
        labels = km.fit_predict(X)
        # silhouette определён только если 1 < k < n_samples
        if k <= 1 or k >= len(X):
            continue
        scores[k] = float(silhouette_score(X, labels))
    return scores


def _best_k(scores: dict[int, float]) -> int:
    if not scores:
        raise ValueError("Не удалось посчитать silhouette ни для одного k.")
    # берём максимум silhouette; при равенстве — меньший k
    return sorted(scores.items(), key=lambda kv: (-kv[1], kv[0]))[0][0]


def _fit_kmeans(X: np.ndarray, k: int, seed: int) -> np.ndarray:
    km = KMeans(n_clusters=k, n_init=50, random_state=seed)
    return km.fit_predict(X)


def _stability_ari(X: np.ndarray, k: int, seeds: tuple[int, ...]) -> tuple[float, float, list[float]]:
    labels_ref = _fit_kmeans(X, k, seeds[0])
    aris: list[float] = []
    for s in seeds[1:]:
        labels = _fit_kmeans(X, k, s)
        aris.append(float(adjusted_rand_score(labels_ref, labels)))
    mean = float(np.mean(aris)) if aris else 1.0
    std = float(np.std(aris)) if aris else 0.0
    return mean, std, aris


def _permute_X_marginals(X: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    # Перемешиваем значения В КАЖДОМ столбце отдельно:
    # маргиналы сохраняются, совместная структура рушится.
    Xp = X.copy()
    for j in range(Xp.shape[1]):
        rng.shuffle(Xp[:, j])
    return Xp


def _perm_test_silhouette(X: np.ndarray, k: int, real_sil: float, B: int, seed: int) -> tuple[np.ndarray, float]:
    rng = np.random.default_rng(seed)
    perm = np.empty(B, dtype=float)
    for b in range(B):
        Xp = _permute_X_marginals(X, rng)
        labels_p = _fit_kmeans(Xp, k, seed + 10_000 + b)
        perm[b] = float(silhouette_score(Xp, labels_p))
    # p-value: доля perm >= real (с псевдосчётом)
    p = float((np.sum(perm >= real_sil) + 1.0) / (B + 1.0))
    return perm, p


def _safe_hist_bins(x: np.ndarray, max_bins: int = 30) -> int:
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    n = len(x)
    if n <= 1:
        return 1
    r = float(np.max(x) - np.min(x))
    # Если диапазон почти нулевой — мало бинов
    if r <= 1e-12:
        return 5
    # sqrt(n) — нормальная дефолтная эвристика
    b = int(max(5, min(max_bins, round(math.sqrt(n)))))
    return b


def _plot_clusters(df: pd.DataFrame, X_raw: np.ndarray, labels: np.ndarray, cfg: Config, fname: str) -> None:
    # Рисуем первые 2 фичи в сыром масштабе, чтобы было читаемо.
    x = df[cfg.features[0]].to_numpy(dtype=float)
    y = df[cfg.features[1]].to_numpy(dtype=float)

    plt.figure(figsize=(10, 7))
    for c in sorted(np.unique(labels)):
        m = labels == c
        plt.scatter(x[m], y[m], s=60, alpha=0.9, label=f"cluster {int(c)}")
    plt.title(f"ECA metric space ({cfg.features[0]},{cfg.features[1]}) with KMeans clusters")
    plt.xlabel(f"Mean {cfg.features[0].replace('_mean','').replace('_',' ')}")
    plt.ylabel(f"Mean {cfg.features[1].replace('_mean','').replace('_',' ')}")
    plt.legend()
    plt.tight_layout()
    plt.savefig(fname, dpi=150)
    plt.close()


def _plot_perm_hist(perm: np.ndarray, real: float, p: float, fname: str) -> None:
    bins = _safe_hist_bins(perm, max_bins=30)

    plt.figure(figsize=(10, 7))
    plt.hist(perm, bins=bins)
    plt.axvline(real, linewidth=2)
    plt.title(f"Permutation test p={p:.4f}")
    plt.xlabel("silhouette score")
    plt.ylabel("count")
    plt.tight_layout()
    plt.savefig(fname, dpi=150)
    plt.close()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default="eca_metrics.csv")
    ap.add_argument("--out-prefix", default="")
    ap.add_argument("--k-min", type=int, default=2)
    ap.add_argument("--k-max", type=int, default=8)
    ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument("--perm-B", type=int, default=500)
    ap.add_argument("--high-entropy-only", action="store_true")
    ap.add_argument("--high-entropy-threshold", type=float, default=0.80)
    ap.add_argument("--features", nargs="+", default=["entropy_mean", "sens_mean"])
    args = ap.parse_args()

    cfg = Config(
        csv_path=args.csv,
        out_prefix=args.out_prefix,
        features=tuple(args.features),
        k_min=args.k_min,
        k_max=args.k_max,
        seed=args.seed,
        stability_seeds=(0, 1, 2, 3, 4, 5, 6, 7, 8, 9),
        perm_B=args.perm_B,
        high_entropy_only=bool(args.high_entropy_only),
        high_entropy_threshold=float(args.high_entropy_threshold),
    )

    df0 = _load_df(cfg.csv_path)
    df = _select(df0, cfg)

    if cfg.high_entropy_only:
        print(f"High-entropy rules: {len(df)}")

    X = _make_X(df, cfg)

    scores = _silhouette_by_k(X, seed=cfg.seed, k_min=cfg.k_min, k_max=cfg.k_max)
    print(f"Silhouette scores: {scores}")

    k_star = _best_k(scores)
    print(f"Selected k* = {k_star}")

    labels = _fit_kmeans(X, k_star, cfg.seed)

    ari_mean, ari_std, ari_list = _stability_ari(X, k_star, cfg.stability_seeds)
    print(f"Stability over seeds (ARI): mean={ari_mean:.3f}, std={ari_std:.3f}")
    if ari_list:
        print(f"ARI list: {ari_list}")

    real_sil = float(silhouette_score(X, labels))
    print(f"Silhouette real={real_sil:.3f}")

    perm, p = _perm_test_silhouette(X, k_star, real_sil, cfg.perm_B, cfg.seed + 999)
    print(f"Permutation p-value={p:.4f}")

    # Сохранения
    out_csv = "eca_clusters.csv" if not cfg.high_entropy_only else "refined_clusters.csv"
    out_scatter = "clusters_HS.png" if not cfg.high_entropy_only else "refined_clusters_HS.png"
    out_perm = "permtest_silhouette.png" if not cfg.high_entropy_only else "refined_permtest.png"

    df_out = df.copy()
    df_out["cluster"] = labels.astype(int)
    df_out.to_csv(out_csv, index=False)

    _plot_clusters(df_out, X, labels, cfg, out_scatter)
    _plot_perm_hist(perm, real_sil, p, out_perm)

    print(f"Saved: {out_scatter}")
    print(f"Saved: {out_csv}")
    print(f"Saved: {out_perm}")


if __name__ == "__main__":
    main()