from __future__ import annotations
from dataclasses import dataclass
from itertools import product
from typing import List, Tuple, Dict, Iterable, Optional, Set

# ----------------------------
# ECA core
# ----------------------------

def rule_from_number(n: int):
    """
    Returns local rule f(a,b,c) for ECA rule n (0..255).
    Bit order: 111,110,101,100,011,010,001,000
    """
    if not (0 <= n <= 255):
        raise ValueError("rule number must be in [0,255]")
    # index: (a,b,c) -> 0..7 where 111 -> 0, 000 -> 7
    def f(a: int, b: int, c: int) -> int:
        idx = (1 - a) * 4 + (1 - b) * 2 + (1 - c)  # maps 111->0, 000->7
        return (n >> (7 - idx)) & 1
    return f

def step_ring(cfg: List[int], f) -> List[int]:
    """One CA step on a ring (periodic boundary)."""
    L = len(cfg)
    out = [0] * L
    for i in range(L):
        a = cfg[(i - 1) % L]
        b = cfg[i]
        c = cfg[(i + 1) % L]
        out[i] = f(a, b, c)
    return out

def shift_ring(cfg: List[int], v: int) -> List[int]:
    """Shift (rotate) ring by v to the right."""
    L = len(cfg)
    v %= L
    if v == 0:
        return cfg[:]
    return cfg[-v:] + cfg[:-v]

def hamming(a: List[int], b: List[int]) -> int:
    return sum(x != y for x, y in zip(a, b))

# ----------------------------
# Domains (periodic fixed points)
# ----------------------------

@dataclass(frozen=True)
class Domain:
    pattern: Tuple[int, ...]  # minimal period word
    name: str = ""

    @property
    def period(self) -> int:
        return len(self.pattern)

    def bit(self, i: int) -> int:
        return self.pattern[i % self.period]

def minimal_period(word: Tuple[int, ...]) -> Tuple[int, ...]:
    """Return minimal repeating block for a periodic word."""
    n = len(word)
    for p in range(1, n + 1):
        if n % p == 0:
            block = word[:p]
            ok = True
            for i in range(n):
                if word[i] != block[i % p]:
                    ok = False
                    break
            if ok:
                return block
    return word

def enumerate_fixed_point_domains(rule_num: int, max_period: int = 8) -> List[Domain]:
    """
    Brute force periodic words of length <= max_period and keep those that are fixed points:
    F(x) = x on the ring of that length.
    Deduplicate by minimal period + rotation.
    """
    f = rule_from_number(rule_num)
    seen: Set[Tuple[int, ...]] = set()
    out: List[Domain] = []

    for n in range(1, max_period + 1):
        for w in product([0, 1], repeat=n):
            cfg = list(w)
            if step_ring(cfg, f) != cfg:
                continue

            block = minimal_period(tuple(cfg))
            # canonicalize up to rotation
            P = len(block)
            rots = [block[i:] + block[:i] for i in range(P)]
            canon = min(rots)
            if canon in seen:
                continue
            seen.add(canon)

            name = f"fixed:{''.join(map(str, canon))}"
            out.append(Domain(pattern=canon, name=name))

    out.sort(key=lambda d: (d.period, d.pattern))
    return out

# ----------------------------
# Two-domain background + defect embedding
# ----------------------------

def background_two_domains(L: int, left: Domain, right: Domain, cut: int) -> List[int]:
    """
    Two-domain background on ring:
    positions [0..cut-1] use left domain, [cut..L-1] use right domain.
    """
    bg = [0] * L
    for i in range(L):
        if i < cut:
            bg[i] = left.bit(i)
        else:
            # phase right domain from cut
            bg[i] = right.bit(i - cut)
    return bg

def embed_defect(bg: List[int], cut: int, defect: Tuple[int, ...]) -> List[int]:
    """
    Replace a segment centered around cut with defect bits.
    Segment starts at cut - len(defect)//2.
    """
    L = len(bg)
    d = len(defect)
    start = (cut - (d // 2)) % L
    cfg = bg[:]
    for j, bit in enumerate(defect):
        cfg[(start + j) % L] = bit
    return cfg

def outside_defect_matches(cfg: List[int], bg: List[int], cut: int, defect_len: int, margin: int) -> bool:
    """
    Locality test: outside [start-margin, start+defect_len+margin) around defect segment,
    cfg must match bg.
    """
    L = len(cfg)
    d = defect_len
    start = (cut - (d // 2)) % L
    # build forbidden set of indices near defect
    bad = set((start + k) % L for k in range(-margin, d + margin))
    for i in range(L):
        if i in bad:
            continue
        if cfg[i] != bg[i]:
            return False
    return True

# ----------------------------
# Particle search (torus approximation)
# ----------------------------

@dataclass(frozen=True)
class ParticleHit:
    left: Domain
    right: Domain
    defect: Tuple[int, ...]
    p: int
    v: int
    mismatch: int

def search_particles_rule69(
    left: Domain,
    right: Domain,
    L: int = 120,
    cut: Optional[int] = None,
    max_defect: int = 8,
    max_p: int = 40,
    margin: int = 4,
    max_mismatch: int = 0,
) -> List[ParticleHit]:
    """
    Search for (p,v) such that after p steps:
      F^p(cfg0) == shift(cfg0, v)
    AND configuration stays localized near the interface (simple locality check).

    Parameters:
      - margin: how many cells around defect are allowed to differ from background
      - max_mismatch: allow small mismatch if you want to catch near-particles
    """
    if cut is None:
        cut = L // 2

    rule_num = 69
    f = rule_from_number(rule_num)
    hits: List[ParticleHit] = []

    bg0 = background_two_domains(L, left, right, cut)

    for d in range(1, max_defect + 1):
        for defect in product([0, 1], repeat=d):
            cfg0 = embed_defect(bg0, cut, defect)

            # enforce locality at t=0 (relative to bg0)
            if not outside_defect_matches(cfg0, bg0, cut, d, margin):
                continue

            cfg = cfg0
            for p in range(1, max_p + 1):
                cfg = step_ring(cfg, f)

                # try all feasible shifts |v|<=p (speed bound for radius 1)
                for v in range(-p, p + 1):
                    target = shift_ring(cfg0, v)
                    mm = hamming(cfg, target)
                    if mm > max_mismatch:
                        continue

                    # locality after p steps: compare to shifted background too
                    bg_shift = shift_ring(bg0, v)
                    if not outside_defect_matches(cfg, bg_shift, cut, d, margin):
                        continue

                    hits.append(ParticleHit(
                        left=left, right=right, defect=tuple(defect),
                        p=p, v=v, mismatch=mm
                    ))

    return hits

# ----------------------------
# Demo runner
# ----------------------------

def main():
    rule_num = 69
    domains = enumerate_fixed_point_domains(rule_num, max_period=8)
    print("Fixed-point periodic domains (Rule 69, period<=8):")
    for d in domains:
        print(f"  period={d.period:2d}  {d.name}")

    # pick two known domains if present:
    # e.g. 01 and 001 (we expect them)
    def pick(pattern_str: str) -> Domain:
        pat = tuple(int(ch) for ch in pattern_str)
        for d in domains:
            if d.pattern == pat:
                return d
        raise RuntimeError(f"domain {pattern_str} not found in period<=8 search")

    left = pick("01")
    right = pick("001")

    print(f"\nSearching particles between {left.name} and {right.name} ...")
    hits = search_particles_rule69(
        left=left,
        right=right,
        L=120,
        max_defect=8,
        max_p=40,
        margin=4,
        max_mismatch=0,
    )

    print(f"Hits: {len(hits)}")
    # show a few best hits (small p then small |v|)
    hits_sorted = sorted(hits, key=lambda h: (h.p, abs(h.v), h.defect))
    for h in hits_sorted[:20]:
        print(f"  p={h.p:2d} v={h.v:3d} speed={h.v/h.p: .3f} defect={''.join(map(str,h.defect))} mismatch={h.mismatch}")

if __name__ == "__main__":
    main()