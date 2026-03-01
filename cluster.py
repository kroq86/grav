#!/usr/bin/env python3
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score, adjusted_rand_score
from sklearn.preprocessing import StandardScaler


def load_features(csv_path="eca_metrics.csv"):
    df = pd.read_csv(csv_path)

    required = ["rule", "entropy_mean", "sens_mean", "lyap_mean"]
    for col in required:
        if col not in df.columns:
            raise ValueError(f"Missing column: {col}")

    return df


def refined_subset(df, entropy_threshold=0.8):
    return df[df["entropy_mean"] > entropy_threshold].copy()


def cluster_analysis(df, k_min=2, k_max=6, seed=0):
    X = df[["entropy_mean", "sens_mean", "lyap_mean"]].to_numpy()
    scaler = StandardScaler()
    X = scaler.fit_transform(X)

    scores = {}
    for k in range(k_min, k_max + 1):
        km = KMeans(n_clusters=k, random_state=seed, n_init=50)
        labels = km.fit_predict(X)
        scores[k] = silhouette_score(X, labels)

    best_k = max(scores, key=scores.get)
    print("Silhouette scores:", scores)
    print("Selected k* =", best_k)

    labels = KMeans(n_clusters=best_k, random_state=seed, n_init=50).fit_predict(X)

    return X, labels, best_k


def stability_test(X, k):
    seeds = [0, 1, 2, 3, 4, 5, 10, 20, 50]
    base = KMeans(n_clusters=k, random_state=seeds[0], n_init=50).fit_predict(X)

    aris = []
    for s in seeds[1:]:
        lab = KMeans(n_clusters=k, random_state=s, n_init=50).fit_predict(X)
        aris.append(adjusted_rand_score(base, lab))

    print("Stability ARI mean =", np.mean(aris))
    return aris


def permutation_test(X, k, B=500):
    rng = np.random.default_rng(123)
    real_labels = KMeans(n_clusters=k, random_state=0, n_init=50).fit_predict(X)
    real_sil = silhouette_score(X, real_labels)

    perm = []
    for _ in range(B):
        Xp = X.copy()
        rng.shuffle(Xp[:, 1])  # ломаем совместную структуру
        lab = KMeans(n_clusters=k, random_state=0, n_init=50).fit_predict(Xp)
        perm.append(silhouette_score(Xp, lab))

    perm = np.array(perm)
    pval = (np.sum(perm >= real_sil) + 1) / (B + 1)

    print("Real silhouette =", real_sil)
    print("Permutation p-value =", pval)

    plt.hist(perm, bins=20)
    plt.axvline(real_sil)
    plt.title(f"Permutation test p={pval:.4f}")
    plt.savefig("refined_permtest.png", dpi=250)
    plt.close()

    return real_sil, pval


if __name__ == "__main__":
    df = load_features()

    df_ref = refined_subset(df, entropy_threshold=0.8)
    print("High-entropy rules:", len(df_ref))

    X, labels, k = cluster_analysis(df_ref)

    stability_test(X, k)

    permutation_test(X, k)

    df_ref["cluster"] = labels
    df_ref.to_csv("refined_clusters.csv", index=False)

    print("Saved refined_clusters.csv and refined_permtest.png")