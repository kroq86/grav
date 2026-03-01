from typing import List, Tuple, Set

# ---- slow rule ----
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

# ---- background + locality ----
def background_two_domains(L: int, left_pat: Tuple[int,...], right_pat: Tuple[int,...], cut: int) -> List[int]:
    cfg = [0]*L
    for i in range(L):
        if i < cut:
            cfg[i] = left_pat[i % len(left_pat)]
        else:
            cfg[i] = right_pat[(i-cut) % len(right_pat)]
    return cfg

def embed_defect(bg: List[int], cut: int, defect: Tuple[int,...]) -> List[int]:
    L = len(bg)
    d = len(defect)
    start = (cut - (d//2)) % L
    cfg = bg[:]
    for j, b in enumerate(defect):
        cfg[(start + j) % L] = b
    return cfg

def defect_indices(L: int, cut: int, d: int, margin: int) -> Set[int]:
    start = (cut - (d//2)) % L
    return { (start + k) % L for k in range(-margin, d + margin) }

def outside_matches(cfg: List[int], bg: List[int], bad: Set[int]) -> bool:
    for i in range(len(cfg)):
        if i in bad:
            continue
        if cfg[i] != bg[i]:
            return False
    return True

# ---- main test ----
def main():
    rule_num = 69
    f = rule_from_number(rule_num)

    L = 240
    cut = 120
    margin = 6

    left_pat = (0,1)
    right_pat = (0,0,1)

    # ⬇️ ВСТАВЬ ИЗ direct_search_fast ВЫВОДА:
    p = 0
    v = 0
    defect_str = "0101"  # заменить!
    defect = tuple(int(c) for c in defect_str)
    d = len(defect)

    bg0 = background_two_domains(L, left_pat, right_pat, cut)
    bad0 = defect_indices(L, cut, d, margin)
    cfg0 = embed_defect(bg0, cut, defect)

    print("t=0 differs-in-zone:",
          any(cfg0[i] != bg0[i] for i in bad0))

    print("t=0 outside matches:",
          outside_matches(cfg0, bg0, bad0))

    cfg = cfg0
    for _ in range(p):
        cfg = step_ring(cfg, f)

    target = shift_ring(cfg0, v)

    print("F^p == shift:", cfg == target)

    bgv = shift_ring(bg0, v)
    badv = { (i + v) % L for i in bad0 }

    print("t=p outside matches shifted bg:",
          outside_matches(cfg, bgv, badv))

if __name__ == "__main__":
    main()