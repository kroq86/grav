from __future__ import annotations
from itertools import product
from typing import Tuple, List, Dict, Optional

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

def bits_from_pattern(pattern: Tuple[int, ...], L: int, offset: int = 0) -> int:
    P = len(pattern)
    x = 0
    for i in range(L):
        b = pattern[(i - offset) % P]
        x |= (b & 1) << i
    return x

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

def int_to_row(x: int, L: int) -> str:
    # MSB left for nicer view
    return "".join("#" if (x >> i) & 1 else "." for i in range(L))

# ------------ search ------------
def find_one_translating_witness(
    rule_num: int,
    left_pat: Tuple[int, ...],
    right_pat: Tuple[int, ...],
    shift: int = -1,
    L: int = 240,
    cut: int = 120,
    max_defect: int = 12,
    margin: int = 6,
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

            # locality at t=1 w.r.t shifted background
            bg1 = rot_r(bg0, L, shift)
            zone1 = rot_r(zone, L, shift)
            outside1 = (~zone1) & m
            if ((x1 ^ bg1) & outside1) != 0:
                continue

            return defect, shift
    return None

def render_spacetime(
    rule_num: int,
    left_pat: Tuple[int, ...],
    right_pat: Tuple[int, ...],
    defect: Tuple[int, ...],
    L: int = 240,
    cut: int = 120,
    steps: int = 80,
) -> List[str]:
    bg0 = background_two_domains_bits(L, left_pat, right_pat, cut)
    x = embed_defect_bits(bg0, L, cut, defect)

    rows = []
    for _ in range(steps):
        rows.append(int_to_row(x, L))
        x = step_eca_bit(x, L, rule_num)
    return rows

def verify_translate(rule_num: int, left_pat: Tuple[int, ...], right_pat: Tuple[int, ...], defect: Tuple[int, ...],
                     shift: int, L: int, cut: int) -> bool:
    bg0 = background_two_domains_bits(L, left_pat, right_pat, cut)
    x0 = embed_defect_bits(bg0, L, cut, defect)
    x1 = step_eca_bit(x0, L, rule_num)
    return x1 == rot_r(x0, L, shift)

def main():
    rule_num = 69
    domains = [
        ("01", (0, 1)),
        ("001", (0, 0, 1)),
        ("00101", (0, 0, 1, 0, 1)),
        ("0010101", (0, 0, 1, 0, 1, 0, 1)),
        ("00100101", (0, 0, 1, 0, 0, 1, 0, 1)),
    ]

    L = 160   # smaller for terminal rendering
    cut = 80
    steps = 60
    max_defect = 12
    margin = 6
    shift = -1

    for i in range(len(domains)):
        for j in range(len(domains)):
            if i == j:
                continue
            nameL, patL = domains[i]
            nameR, patR = domains[j]

            w = find_one_translating_witness(rule_num, patL, patR, shift=shift, L=L, cut=cut,
                                             max_defect=max_defect, margin=margin)
            if not w:
                continue

            defect, s = w
            ok = verify_translate(rule_num, patL, patR, defect, s, L, cut)

            print(f"\nPair {nameL} | {nameR}  shift={s:+d}  defect={''.join(map(str,defect))}  VERIFY={ok}")
            rows = render_spacetime(rule_num, patL, patR, defect, L=L, cut=cut, steps=steps)
            # print a compact view (first 30 rows)
            for r in rows[:30]:
                print(r)

            # stop after first pair to avoid spamming terminal
            return

if __name__ == "__main__":
    main()