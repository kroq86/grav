# find_translating_fronts.py
# Fast bitwise search for translating fronts in ECA Rule 69:
#   F(x) = shift(x, s) with s in {-1, +1}
# with locality relative to a two-domain background.
#
# Usage:
#   python3 find_translating_fronts.py
#
# Notes:
# - This searches "front" (translating solutions) rather than particles with p>1.
# - It can iterate over all pairs of periodic fixed-point domains (provided below),
#   or you can restrict to a single pair.

from __future__ import annotations
from itertools import product
from typing import Dict, Tuple, List, Optional

# ----------------------------
# Bit helpers (ring of length L)
# ----------------------------

def mask_L(L: int) -> int:
    return (1 << L) - 1

def rot_r(x: int, L: int, v: int) -> int:
    """Rotate right by v (v can be negative)."""
    v %= L
    m = mask_L(L)
    if v == 0:
        return x & m
    return ((x >> v) | ((x << (L - v)) & m)) & m

# ----------------------------
# ECA step (bitwise)
# ----------------------------

def step_eca_bit(x: int, L: int, rule_num: int) -> int:
    """
    One ECA step on a ring using bit operations.
    Bit i corresponds to cell i.
    Neighborhood (a,b,c) = (left, center, right) where:
      left  = cell i-1
      right = cell i+1
    Wolfram rule bits order: 111,110,101,100,011,010,001,000.
    """
    m = mask_L(L)
    x &= m

    # left neighbor of i is i-1 -> shift right by 1 with wrap
    left = (x >> 1) | ((x & 1) << (L - 1))
    # right neighbor of i is i+1 -> shift left by 1 with wrap
    right = ((x << 1) & m) | (x >> (L - 1))
    center = x

    nxt = 0
    for a in (0, 1):
        La = left if a else (~left)
        for b in (0, 1):
            Cb = center if b else (~center)
            for c in (0, 1):
                Rc = right if c else (~right)
                pat_mask = La & Cb & Rc & m
                idx = (a << 2) | (b << 1) | c  # 000..111
                out = (rule_num >> (7 - idx)) & 1
                if out:
                    nxt |= pat_mask
    return nxt & m

# ----------------------------
# Domains on a ring
# ----------------------------

def bits_from_pattern(pattern: Tuple[int, ...], L: int, offset: int = 0) -> int:
    """Periodic pattern to ring bits of length L, with phase offset."""
    P = len(pattern)
    x = 0
    for i in range(L):
        b = pattern[(i - offset) % P]
        x |= (b & 1) << i
    return x

def is_global_periodic(x: int, L: int, pat: Tuple[int, ...]) -> bool:
    """x equals some rotation of pure domain pattern on ring."""
    base = bits_from_pattern(pat, L, offset=0)
    P = len(pat)
    for phase in range(P):
        if x == rot_r(base, L, phase):
            return True
    return False

# ----------------------------
# Two-domain background + defect zone
# ----------------------------

def background_two_domains_bits(L: int, left_pat: Tuple[int, ...], right_pat: Tuple[int, ...], cut: int) -> int:
    """Left domain on [0..cut-1], right domain on [cut..L-1]."""
    x = 0
    PL = len(left_pat)
    PR = len(right_pat)
    for i in range(L):
        if i < cut:
            b = left_pat[i % PL]
        else:
            b = right_pat[(i - cut) % PR]
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

# ----------------------------
# Translating front search
# ----------------------------

def find_translating_fronts(
    rule_num: int,
    left_pat: Tuple[int, ...],
    right_pat: Tuple[int, ...],
    shifts: Tuple[int, ...] = (-1, +1),
    L: int = 240,
    cut: int = 120,
    max_defect: int = 12,
    margin: int = 6,
    max_witness_per_shift: int = 5,
) -> Dict[int, List[Tuple[int, ...]]]:
    """
    Returns dict: shift s -> list of witness defects (tuples) such that:
      x1 = F(x0) == rot_r(x0, s)
    with locality relative to bg0 and shifted bg0.

    Locality:
      - at t=0: outside defect-zone cfg matches bg0
      - at t=1: outside shifted defect-zone cfg matches shifted bg0
    """
    m = mask_L(L)
    bg0 = background_two_domains_bits(L, left_pat, right_pat, cut)
    out: Dict[int, List[Tuple[int, ...]]] = {s: [] for s in shifts}

    for d in range(1, max_defect + 1):
        zone = build_zone_mask(L, defect_zone_indices(L, cut, d, margin)) & m
        outside = (~zone) & m

        for defect in product([0, 1], repeat=d):
            defect = tuple(defect)
            x0 = embed_defect_bits(bg0, L, cut, defect)

            # must differ from bg inside zone (real defect)
            if ((x0 ^ bg0) & zone) == 0:
                continue

            # must match bg outside zone at t=0
            if ((x0 ^ bg0) & outside) != 0:
                continue

            # reject trivial global pure domains
            if is_global_periodic(x0, L, left_pat) or is_global_periodic(x0, L, right_pat):
                continue

            x1 = step_eca_bit(x0, L, rule_num)

            for s in shifts:
                if len(out[s]) >= max_witness_per_shift:
                    continue

                # exact translate per step
                if x1 != rot_r(x0, L, s):
                    continue

                # locality at t=1 relative to shifted background
                bg1 = rot_r(bg0, L, s)
                zone1 = rot_r(zone, L, s)
                outside1 = (~zone1) & m
                if ((x1 ^ bg1) & outside1) != 0:
                    continue

                out[s].append(defect)

        # early stop if we have enough witnesses for all shifts
        if all(len(out[s]) >= max_witness_per_shift for s in shifts):
            break

    return out

# ----------------------------
# Run: all pairs of known fixed-point domains for Rule 69
# (from your earlier output)
# ----------------------------

def main():
    rule_num = 69

    # Fixed-point periodic domains (period<=8) you listed:
    domains = [
        ("01", (0, 1)),
        ("001", (0, 0, 1)),
        ("00101", (0, 0, 1, 0, 1)),
        ("0010101", (0, 0, 1, 0, 1, 0, 1)),
        ("00100101", (0, 0, 1, 0, 0, 1, 0, 1)),
    ]

    # Parameters
    L = 240
    cut = 120
    max_defect = 12
    margin = 6

    print(f"Rule {rule_num} translating fronts search")
    print(f"L={L} cut={cut} max_defect={max_defect} margin={margin}\n")

    for i in range(len(domains)):
        for j in range(len(domains)):
            if i == j:
                continue
            nameL, patL = domains[i]
            nameR, patR = domains[j]

            res = find_translating_fronts(
                rule_num,
                patL,
                patR,
                shifts=(-1, +1),
                L=L,
                cut=cut,
                max_defect=max_defect,
                margin=margin,
                max_witness_per_shift=3,
            )

            if not res[-1] and not res[+1]:
                continue

            print(f"Pair {nameL} | {nameR}:")
            for s in (-1, +1):
                if res[s]:
                    for d in res[s]:
                        dstr = "".join(map(str, d))
                        print(f"  shift={s:+d}  witness_defect={dstr}")
            print()

if __name__ == "__main__":
    main()