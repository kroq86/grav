from __future__ import annotations
from dataclasses import dataclass
from itertools import product
from typing import List, Tuple, Dict, Set

# ---------- ECA ----------
def rule_from_number(n: int):
    def f(a: int, b: int, c: int) -> int:
        idx = (1 - a) * 4 + (1 - b) * 2 + (1 - c)
        return (n >> (7 - idx)) & 1
    return f

def step_ring(cfg: List[int], f) -> List[int]:
    L = len(cfg)
    out = [0] * L
    for i in range(L):
        out[i] = f(cfg[(i-1) % L], cfg[i], cfg[(i+1) % L])
    return out

def shift_ring(cfg: List[int], v: int) -> List[int]:
    L = len(cfg)
    v %= L
    return cfg[-v:] + cfg[:-v] if v else cfg[:]

def hamming(a: List[int], b: List[int]) -> int:
    return sum(x != y for x, y in zip(a, b))

# ---------- Domains ----------
@dataclass(frozen=True)
class Domain:
    pattern: Tuple[int, ...]
    name: str

    @property
    def P(self) -> int:
        return len(self.pattern)

    def bit(self, i: int) -> int:
        return self.pattern[i % self.P]

def is_global_domain(cfg: List[int], d: Domain) -> bool:
    """cfg equals some phase rotation of domain d (exact periodic)."""
    L = len(cfg)
    P = d.P
    if L % P != 0:
        return False
    for phase in range(P):
        ok = True
        for i in range(L):
            if cfg[i] != d.bit(i - phase):
                ok = False
                break
        if ok:
            return True
    return False

# ---------- Background + locality ----------
def background_two_domains(L: int, DL: Domain, DR: Domain, cut: int) -> List[int]:
    bg = [0]*L
    for i in range(L):
        if i < cut:
            bg[i] = DL.bit(i)
        else:
            bg[i] = DR.bit(i - cut)
    return bg

def defect_indices(L: int, cut: int, d: int, margin: int) -> Set[int]:
    start = (cut - (d // 2)) % L
    return set((start + k) % L for k in range(-margin, d + margin))

def outside_matches(cfg: List[int], bg: List[int], bad: Set[int]) -> bool:
    for i in range(len(cfg)):
        if i in bad:
            continue
        if cfg[i] != bg[i]:
            return False
    return True

def embed_defect(bg: List[int], cut: int, defect: Tuple[int, ...]) -> List[int]:
    L = len(bg)
    d = len(defect)
    start = (cut - (d // 2)) % L
    cfg = bg[:]
    for j, b in enumerate(defect):
        cfg[(start + j) % L] = b
    return cfg

# ---------- Search ----------
@dataclass(frozen=True)
class Hit:
    p: int
    v: int
    defect: Tuple[int, ...]

def search_direct(
    rule_num: int,
    DL: Domain,
    DR: Domain,
    L: int = 240,
    cut: int = 120,
    max_defect: int = 10,
    max_p: int = 60,
    margin: int = 6,
) -> Dict[Tuple[int,int], List[Tuple[int,...]]]:

    f = rule_from_number(rule_num)
    bg0 = background_two_domains(L, DL, DR, cut)

    results: Dict[Tuple[int,int], List[Tuple[int,...]]] = {}

    for d in range(1, max_defect+1):
        bad0 = defect_indices(L, cut, d, margin)

        for defect in product([0,1], repeat=d):
            cfg0 = embed_defect(bg0, cut, defect)

            # must be a real defect (differs from bg in allowed zone)
            if not any(cfg0[i] != bg0[i] for i in bad0):
                continue

            # reject if cfg0 became a global pure domain (trivial)
            if is_global_domain(cfg0, DL) or is_global_domain(cfg0, DR):
                continue

            # enforce locality at t=0
            if not outside_matches(cfg0, bg0, bad0):
                continue

            cfg = cfg0
            for p in range(2, max_p+1):  # p>=2
                cfg = step_ring(cfg, f)

                # speed bound radius 1: |v| <= p
                for v in range(-p, p+1):
                    if v == 0:
                        continue

                    if cfg != shift_ring(cfg0, v):
                        continue

                    # locality after p steps relative to shifted background
                    bgv = shift_ring(bg0, v)
                    badv = set((i + v) % L for i in bad0)  # defect zone shifts too
                    if not outside_matches(cfg, bgv, badv):
                        continue

                    results.setdefault((p,v), []).append(tuple(defect))

    return results

if __name__ == "__main__":
    rule_num = 69
    DL = Domain((0,1), "01")
    DR = Domain((0,0,1), "001")

    res = search_direct(
        rule_num, DL, DR,
        L=240, cut=120,
        max_defect=10,
        max_p=60,
        margin=6
    )

    print(f"Found (p,v) pairs: {len(res)}")
    for (p,v) in sorted(res.keys(), key=lambda x: (x[0], abs(x[1]), x[1]))[:30]:
        print(f"  p={p:2d} v={v:3d} speed={v/p: .3f}  examples={min(5, len(res[(p,v)]))}")