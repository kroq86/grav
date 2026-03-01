from __future__ import annotations
from typing import List, Tuple
from itertools import product

# ---- slow ----
def rule_from_number(n: int):
    def f(a: int, b: int, c: int) -> int:
        idx = (1 - a) * 4 + (1 - b) * 2 + (1 - c)
        return (n >> (7 - idx)) & 1
    return f

def step_ring(cfg: List[int], f) -> List[int]:
    L = len(cfg)
    out = [0]*L
    for i in range(L):
        out[i] = f(cfg[(i-1)%L], cfg[i], cfg[(i+1)%L])
    return out

def shift_ring(cfg: List[int], v: int) -> List[int]:
    L = len(cfg)
    v %= L
    return cfg[-v:] + cfg[:-v] if v else cfg[:]

# ---- fast (copy exact from direct_search_fast.py!) ----
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
    # Keep the same cell indexing convention as nonlinear.py/make_rule:
    # bit i is cell i, left neighbor is i-1, right neighbor is i+1.
    left = ((x << 1) & m) | (x >> (L - 1))
    right = (x >> 1) | ((x & 1) << (L - 1))
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
                out = (rule_num >> idx) & 1
                if out:
                    nxt |= pat_mask
    return nxt & m

def int_to_list(x: int, L: int) -> List[int]:
    return [(x >> i) & 1 for i in range(L)]

def list_to_int(cfg: List[int]) -> int:
    x = 0
    for i, b in enumerate(cfg):
        if b:
            x |= 1 << i
    return x

def background_two_domains_list(L: int, left_pat: Tuple[int,...], right_pat: Tuple[int,...], cut: int) -> List[int]:
    cfg = [0]*L
    for i in range(L):
        cfg[i] = left_pat[i % len(left_pat)] if i < cut else right_pat[(i-cut) % len(right_pat)]
    return cfg

def embed_defect_list(bg: List[int], cut: int, defect: Tuple[int,...]) -> List[int]:
    L = len(bg)
    d = len(defect)
    start = (cut - (d//2)) % L
    cfg = bg[:]
    for j, b in enumerate(defect):
        cfg[(start + j) % L] = b
    return cfg

def main():
    rule_num = 69
    L = 240
    cut = 120
    left_pat = (0,1)
    right_pat = (0,0,1)

    # 👇 ВСТАВЬ СЮДА значения из direct_search_fast.py вывода:
    p = 0
    v = 0
    defect_str = "0101"  # пример
    defect = tuple(int(ch) for ch in defect_str)

    bg = background_two_domains_list(L, left_pat, right_pat, cut)
    x0_list = embed_defect_list(bg, cut, defect)
    x0_int = list_to_int(x0_list)

    f = rule_from_number(rule_num)

    # evolve slow
    xs = x0_list
    for _ in range(p):
        xs = step_ring(xs, f)

    # evolve fast
    xf = x0_int
    for _ in range(p):
        xf = step_eca_bit(xf, L, rule_num)
    xf_list = int_to_list(xf, L)

    print("slow==fast after p steps:", xs == xf_list)

    # check equality to shift
    ok_slow = (xs == shift_ring(x0_list, v))
    ok_fast = (xf == rot_r(x0_int, L, v))
    print("slow verifies F^p=shift:", ok_slow)
    print("fast verifies F^p=shift:", ok_fast)

if __name__ == "__main__":
    main()
