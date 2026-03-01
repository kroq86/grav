from __future__ import annotations
from itertools import product
from typing import Tuple, List, Dict, Optional
import os, json
import matplotlib.pyplot as plt

# ------------ bit CA core ------------
def mask_L(L: int) -> int:
    return (1 << L) - 1

def rot_r(x: int, L: int, v: int) -> int:
    v %= L
    m = mask_L(L)
    if v == 0:
        return x & m
    return ((x >> v) | ((x << (L - v)) & m)) & m

def step_eca_bit(x: int, L: int, rule_num: int) -> int:
    m = mask_L(L)
    x &= m
    left = (x >> 1) | ((x & 1) << (L - 1))            # i-1
    right = ((x << 1) & m) | (x >> (L - 1))           # i+1
    center = x

    nxt = 0
    for a in (0, 1):
        La = left if a else (~left)
        for b in (0, 1):
            Cb = center if b else (~center)
            for c in (0, 1):
                Rc = right if c else (~right)
                pat = La & Cb & Rc & m
                idx = (a << 2) | (b << 1) | c
                out = (rule_num >> (7 - idx)) & 1
                if out:
                    nxt |= pat
    return nxt & m

def background_two_domains_bits(L: int, left_pat: Tuple[int, ...], right_pat: Tuple[int, ...], cut: int) -> int:
    x = 0
    PL = len(left_pat)
    PR = len(right_pat)
    for i in range(L):
        b = left_pat[i % PL] if i < cut else right_pat[(i - cut) % PR]
        x |= (b & 1) << i
    return x

def defect_zone_indices(L: int, cut: int, d: int, margin: int) -> List[int]:
    start = (cut - (d // 2)) % L
    return [(start + k) % L for k in range(-margin, d + margin)]

def build_zone_mask(L: int, idxs: List[int]) -> int:
    z = 0
    for i in idxs:
        z |= 1 << i
    return z

def embed_defect_bits(bg: int, L: int, cut: int, defect: Tuple[int, ...]) -> int:
    d = len(defect)
    start = (cut - (d // 2)) % L
    x = bg
    for j, b in enumerate(defect):
        i = (start + j) % L
        if b:
            x |= (1 << i)
        else:
            x &= ~(1 << i)
    return x & mask_L(L)

# ------------ search witness ------------
def find_one_translating_witness(
    rule_num: int,
    left_pat: Tuple[int, ...],
    right_pat: Tuple[int, ...],
    shift: int,
    L: int,
    cut: int,
    max_defect: int,
    margin: int,
) -> Optional[Tuple[Tuple[int, ...], int]]:
    m = mask_L(L)
    bg0 = background_two_domains_bits(L, left_pat, right_pat, cut)

    for d in range(1, max_defect + 1):
        zone = build_zone_mask(L, defect_zone_indices(L, cut, d, margin)) & m
        outside = (~zone) & m

        for defect in product([0, 1], repeat=d):
            defect = tuple(defect)
            x0 = embed_defect_bits(bg0, L, cut, defect)

            if ((x0 ^ bg0) & zone) == 0:
                continue
            if ((x0 ^ bg0) & outside) != 0:
                continue

            x1 = step_eca_bit(x0, L, rule_num)
            if x1 != rot_r(x0, L, shift):
                continue

            bg1 = rot_r(bg0, L, shift)
            zone1 = rot_r(zone, L, shift)
            outside1 = (~zone1) & m
            if ((x1 ^ bg1) & outside1) != 0:
                continue

            return defect, shift
    return None

def spacetime_matrix(
    rule_num: int,
    left_pat: Tuple[int, ...],
    right_pat: Tuple[int, ...],
    defect: Tuple[int, ...],
    L: int,
    cut: int,
    steps: int,
) -> List[List[int]]:
    bg0 = background_two_domains_bits(L, left_pat, right_pat, cut)
    x = embed_defect_bits(bg0, L, cut, defect)
    mat: List[List[int]] = []
    for _ in range(steps):
        row = [(x >> i) & 1 for i in range(L)]
        mat.append(row)
        x = step_eca_bit(x, L, rule_num)
    return mat

def save_png(mat: List[List[int]], path: str):
    # mat[t][i] with i increasing to the right
    # show i from 0..L-1 left->right (imshow expects x-axis increasing)
    import numpy as np
    arr = np.array(mat, dtype=int)
    plt.figure(figsize=(12, 6))
    plt.imshow(arr, interpolation="nearest", aspect="auto", cmap="gray_r")
    plt.xlabel("cell index")
    plt.ylabel("time")
    plt.tight_layout()
    plt.savefig(path, dpi=200)
    plt.close()

def main():
    rule_num = 69
    domains = [
        ("01", (0, 1)),
        ("001", (0, 0, 1)),
        ("00101", (0, 0, 1, 0, 1)),
        ("0010101", (0, 0, 1, 0, 1, 0, 1)),
        ("00100101", (0, 0, 1, 0, 0, 1, 0, 1)),
    ]

    outdir = "fronts_out"
    os.makedirs(outdir, exist_ok=True)

    L = 240
    cut = 120
    steps = 160
    max_defect = 12
    margin = 6
    shift = -1

    certs = []

    for i in range(len(domains)):
        for j in range(len(domains)):
            if i == j:
                continue

            nameL, patL = domains[i]
            nameR, patR = domains[j]

            w = find_one_translating_witness(rule_num, patL, patR, shift, L, cut, max_defect, margin)
            if not w:
                continue

            defect, s = w
            dstr = "".join(map(str, defect))

            mat = spacetime_matrix(rule_num, patL, patR, defect, L, cut, steps)
            png_name = f"front_rule{rule_num}_{nameL}_to_{nameR}_shift{s:+d}_def{dstr}.png"
            png_path = os.path.join(outdir, png_name)
            save_png(mat, png_path)

            cert = {
                "rule": rule_num,
                "left_domain": nameL,
                "right_domain": nameR,
                "shift_per_step": s,
                "defect_bits": dstr,
                "L": L,
                "cut": cut,
                "steps_rendered": steps,
                "max_defect_searched": max_defect,
                "margin": margin,
            }
            certs.append(cert)

    with open(os.path.join(outdir, "certificates.json"), "w", encoding="utf-8") as f:
        json.dump(certs, f, ensure_ascii=False, indent=2)

    print(f"Saved {len(certs)} PNGs + certificates.json into {outdir}/")

if __name__ == "__main__":
    main()