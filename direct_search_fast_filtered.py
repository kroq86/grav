from __future__ import annotations
from itertools import product
from typing import Dict, Tuple, List, Set

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
    # left neighbor i-1
    left = (x >> 1) | ((x & 1) << (L - 1))
    # right neighbor i+1
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
                idx = (a << 2) | (b << 1) | c
                out = (rule_num >> (7 - idx)) & 1
                if out:
                    nxt |= pat_mask
    return nxt & m

def bits_from_pattern(pattern: Tuple[int, ...], L: int, offset: int = 0) -> int:
    P = len(pattern)
    x = 0
    for i in range(L):
        b = pattern[(i - offset) % P]
        x |= (b & 1) << i
    return x

def is_global_periodic(x: int, L: int, pat: Tuple[int, ...]) -> bool:
    base = bits_from_pattern(pat, L, offset=0)
    P = len(pat)
    for phase in range(P):
        if x == rot_r(base, L, phase):
            return True
    return False

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
    return [ (start + k) % L for k in range(-margin, d + margin) ]

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
            x |= 1 << i
        else:
            x &= ~(1 << i)
    return x & mask_L(L)

def first_shift_if_translate_per_step(x0: int, L: int, rule_num: int) -> Tuple[bool, int]:
    """Check if F(x0) is a pure shift of x0. If yes, return (True, v)."""
    x1 = step_eca_bit(x0, L, rule_num)
    # For radius 1, only shifts in {-1,0,+1} are plausible per step for localized stuff,
    # but translating solutions can still be checked in wider range if you want.
    for v in (-1, +1):
        if x1 == rot_r(x0, L, v):
            return True, v
    return False, 0

def search_filtered(
    rule_num: int,
    left_pat: Tuple[int, ...],
    right_pat: Tuple[int, ...],
    L: int = 240,
    cut: int = 120,
    max_defect: int = 12,
    max_p: int = 120,
    margin: int = 6,
) -> Dict[Tuple[int, int], Tuple[int, ...]]:
    m = mask_L(L)
    bg0 = background_two_domains_bits(L, left_pat, right_pat, cut)
    results: Dict[Tuple[int, int], Tuple[int, ...]] = {}

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
            if is_global_periodic(x0, L, left_pat) or is_global_periodic(x0, L, right_pat):
                continue

            # FILTER: if translates by a shift each single step, treat as "front" and skip
            ok1, v1 = first_shift_if_translate_per_step(x0, L, rule_num)
            if ok1:
                continue

            x = x0
            for p in range(2, max_p + 1):
                x = step_eca_bit(x, L, rule_num)

                for v in range(-p, p + 1):
                    if v == 0:
                        continue
                    if x != rot_r(x0, L, v):
                        continue

                    bgv = rot_r(bg0, L, v)
                    zonev = rot_r(zone, L, v)
                    outsidev = (~zonev) & m
                    if ((x ^ bgv) & outsidev) != 0:
                        continue

                    results.setdefault((p, v), defect)

    return results

if __name__ == "__main__":
    rule_num = 69
    left_pat = (0, 1)
    right_pat = (0, 0, 1)

    res = search_filtered(
        rule_num, left_pat, right_pat,
        L=240, cut=120,
        max_defect=12,
        max_p=120,
        margin=6
    )

    print(f"Filtered (p,v) pairs: {len(res)}")
    for (p, v) in sorted(res.keys(), key=lambda x: (x[0], abs(x[1]), x[1]))[:60]:
        defect = res[(p, v)]
        print(f"  p={p:3d} v={v:4d} speed={v/p: .6f} defect={''.join(map(str, defect))}")