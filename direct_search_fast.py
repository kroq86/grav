from __future__ import annotations
from itertools import product
from typing import Dict, Tuple, List

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
# Build periodic domains on a ring
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
    """
    Left domain on [0..cut-1], right domain on [cut..L-1].
    """
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
            x |= (1 << i)
        else:
            x &= ~(1 << i)
    return x & mask_L(L)

# ----------------------------
# Fast direct search
# ----------------------------

def search_direct_fast(
    rule_num: int,
    left_pat: Tuple[int, ...],
    right_pat: Tuple[int, ...],
    L: int = 240,
    cut: int = 120,
    max_defect: int = 12,
    max_p: int = 120,
    margin: int = 8,
) -> Dict[Tuple[int, int], Tuple[int, ...]]:
    """
    Returns dict (p,v) -> witness defect (tuple of 0/1).
    Conditions:
      - p >= 2
      - v != 0
      - |v| <= p (radius 1 light-cone bound)
      - exact equality: F^p(x0) == rot_r(x0, v)
      - locality: outside defect zone matches background (t=0 and t=p w.r.t shifted bg)
      - rejects trivial pure-domain configs
    """
    m = mask_L(L)
    bg0 = background_two_domains_bits(L, left_pat, right_pat, cut)

    results: Dict[Tuple[int, int], Tuple[int, ...]] = {}

    for d in range(1, max_defect + 1):
        zone_idxs = defect_zone_indices(L, cut, d, margin)
        zone = build_zone_mask(L, zone_idxs) & m
        outside = (~zone) & m

        for defect in product([0, 1], repeat=d):
            defect = tuple(defect)
            x0 = embed_defect_bits(bg0, L, cut, defect)

            # must differ from bg inside zone
            if ((x0 ^ bg0) & zone) == 0:
                continue

            # must match bg outside zone at t=0
            if ((x0 ^ bg0) & outside) != 0:
                continue

            # reject trivial global pure domains
            if is_global_periodic(x0, L, left_pat) or is_global_periodic(x0, L, right_pat):
                continue

            x = x0
            for p in range(2, max_p + 1):
                x = step_eca_bit(x, L, rule_num)

                # |v|<=p, v!=0
                for v in range(-p, p + 1):
                    if v == 0:
                        continue

                    if x != rot_r(x0, L, v):
                        continue

                    # locality after p relative to shifted background
                    bgv = rot_r(bg0, L, v)
                    zonev = rot_r(zone, L, v)
                    outsidev = (~zonev) & m

                    if ((x ^ bgv) & outsidev) != 0:
                        continue

                    # store first witness for (p,v)
                    results.setdefault((p, v), defect)

    return results

# ----------------------------
# CLI demo
# ----------------------------

if __name__ == "__main__":
    rule_num = 69
    left_pat = (0, 1)        # 01
    right_pat = (0, 0, 1)    # 001

    res = search_direct_fast(
        rule_num,
        left_pat,
        right_pat,
        L=240,
        cut=120,
        max_defect=12,
        max_p=120,
        margin=6,
    )

    print(f"Found (p,v) pairs: {len(res)}")
    for (p, v) in sorted(res.keys(), key=lambda x: (x[0], abs(x[1]), x[1]))[:60]:
        defect = res[(p, v)]
        dstr = "".join(map(str, defect))
        print(f"  p={p:3d} v={v:4d} speed={v/p: .6f}  defect={dstr}")