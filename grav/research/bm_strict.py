from __future__ import annotations
from dataclasses import dataclass
from itertools import product
from typing import Tuple, Dict, List, Set, Optional

# ---------- ECA ----------
def rule_from_number(n: int):
    def f(a: int, b: int, c: int) -> int:
        idx = (1 - a) * 4 + (1 - b) * 2 + (1 - c)  # 111->0 ... 000->7
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

# ---------- Domains ----------
@dataclass(frozen=True)
class Domain:
    pattern: Tuple[int, ...]
    name: str = ""

    @property
    def P(self) -> int:
        return len(self.pattern)

    def bit(self, i: int) -> int:
        return self.pattern[i % self.P]

# ---------- Boundary machine ----------
@dataclass(frozen=True)
class State:
    phase_L: int
    phase_R: int
    window: Tuple[int, ...]  # length W

@dataclass(frozen=True)
class Edge:
    to: State
    shift: int  # -1,0,+1

class BoundaryMachine:
    def __init__(self, rule_num: int, DL: Domain, DR: Domain, W: int):
        if W % 2 == 0:
            raise ValueError("Use odd W (e.g. 7,9,11) for a centered boundary.")
        self.f = rule_from_number(rule_num)
        self.rule_num = rule_num
        self.DL = DL
        self.DR = DR
        self.W = W
        self.states: Set[State] = set()
        self.graph: Dict[State, List[Edge]] = {}

    def background_window(self, phase_L: int, phase_R: int) -> Tuple[int, ...]:
        W = self.W
        mid = W // 2
        bg = []
        for i in range(W):
            if i < mid:
                bg.append(self.DL.bit((i - mid) + phase_L))
            else:
                bg.append(self.DR.bit((i - mid) + phase_R))
        return tuple(bg)

    def initial_states(self) -> List[State]:
        W = self.W
        bg00 = self.background_window(0, 0)
        inits = []
        for w in product([0, 1], repeat=W):
            w = tuple(w)
            if w == bg00:
                continue
            inits.append(State(0, 0, w))
        return inits

    @staticmethod
    def _hamming_tuple(a: Tuple[int, ...], b: Tuple[int, ...]) -> int:
        return sum(x != y for x, y in zip(a, b))

    def _best_shift(self, new_window: Tuple[int, ...], phase_L: int, phase_R: int) -> int:
        """
        Choose shift in {-1,0,+1} that best matches the ideal background window
        after applying that shift (and updating phases accordingly).
        """
        candidates = []
        for sh in (-1, 0, +1):
            pl = (phase_L + sh) % self.DL.P
            pr = (phase_R + sh) % self.DR.P
            bg = self.background_window(pl, pr)
            score = self._hamming_tuple(new_window, bg)
            candidates.append((score, abs(sh), sh))  # tie-break: prefer 0 shift
        candidates.sort()
        return candidates[0][2]

    def step_state(self, s: State) -> Edge:
        W = self.W

        # reconstruct extended config length W+2 with domains outside
        ext = []
        for i in range(-1, W + 1):
            if 0 <= i < W:
                ext.append(s.window[i])
            else:
                if i < 0:
                    ext.append(self.DL.bit(i + s.phase_L))
                else:
                    ext.append(self.DR.bit((i - W) + s.phase_R))

        # apply local rule to get new window
        new = []
        for i in range(1, W + 1):
            new.append(self.f(ext[i - 1], ext[i], ext[i + 1]))
        new = tuple(new)

        # STRICT shift: argmin mismatch with ideal background
        sh = self._best_shift(new, s.phase_L, s.phase_R)

        new_phase_L = (s.phase_L + sh) % self.DL.P
        new_phase_R = (s.phase_R + sh) % self.DR.P

        return Edge(State(new_phase_L, new_phase_R, new), sh)

    def build(self, max_states: int = 300000):
        stack = self.initial_states()
        while stack and len(self.states) < max_states:
            s = stack.pop()
            if s in self.states:
                continue
            self.states.add(s)
            e = self.step_state(s)
            self.graph.setdefault(s, []).append(e)
            if e.to not in self.states:
                stack.append(e.to)

    def find_cycles(self, max_len: int = 80) -> List[Tuple[List[State], List[int]]]:
        """
        Return cycles as (states, shifts), where states includes start repeated at end.
        """
        cycles = []
        for start in self.states:
            stack = [(start, [start], [])]
            while stack:
                cur, path, shifts = stack.pop()
                if len(path) > max_len:
                    continue
                for e in self.graph.get(cur, []):
                    if e.to == start and len(path) >= 2:
                        cycles.append((path + [start], shifts + [e.shift]))
                    elif e.to not in path:
                        stack.append((e.to, path + [e.to], shifts + [e.shift]))
        return cycles

# ---------- Witness build + verify ----------
def build_witness_on_ring(DL: Domain, DR: Domain, W: int, state: State, L: int = 240, cut: Optional[int] = None) -> List[int]:
    if cut is None:
        cut = L // 2
    mid = W // 2
    cfg = [0] * L
    for i in range(L):
        if i < cut:
            cfg[i] = DL.bit(i + state.phase_L)
        else:
            cfg[i] = DR.bit((i - cut) + state.phase_R)
    start = (cut - mid) % L
    for j, b in enumerate(state.window):
        cfg[(start + j) % L] = b
    return cfg

def verify_shift(rule_num: int, cfg0: List[int], p: int, v: int) -> bool:
    f = rule_from_number(rule_num)
    cfg = cfg0[:]
    for _ in range(p):
        cfg = step_ring(cfg, f)
    return cfg == shift_ring(cfg0, v)

# ---------- Run ----------
if __name__ == "__main__":
    rule_num = 69
    DL = Domain((0, 1), "01")
    DR = Domain((0, 0, 1), "001")

    W = 9  # start at 9; try 7/9/11
    bm = BoundaryMachine(rule_num, DL, DR, W)
    bm.build()

    cycles = bm.find_cycles(max_len=60)
    print(f"States: {len(bm.states)}  cycles found (raw): {len(cycles)}")

    verified = []
    for states, shifts in cycles:
        p = len(shifts)
        v = sum(shifts)
        if v == 0:
            continue
        s0 = states[0]
        cfg0 = build_witness_on_ring(DL, DR, W, s0, L=300)
        if verify_shift(rule_num, cfg0, p, v):
            verified.append((p, v, s0))

    verified.sort(key=lambda x: (x[0], abs(x[1]), x[1]))
    print("Verified moving cycles (p,v):")
    for p, v, _ in verified[:20]:
        print(f"  p={p:2d} v={v:3d} speed={v/p: .3f}")

    if not verified:
        print("No verified moving cycles at this W. Try increasing W to 11 or 13.")