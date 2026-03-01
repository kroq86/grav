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
    shift: int

class BoundaryMachine:
    def __init__(self, rule_num: int, DL: Domain, DR: Domain, W: int):
        self.f = rule_from_number(rule_num)
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
        inits = []
        bg00 = self.background_window(0, 0)
        for w in product([0, 1], repeat=W):
            if tuple(w) == bg00:
                continue
            inits.append(State(0, 0, tuple(w)))
        return inits

    def step_state(self, s: State) -> Edge:
        W = self.W
        mid = W // 2

        # reconstruct extended config of length W+2
        ext = []
        for i in range(-1, W+1):
            if 0 <= i < W:
                ext.append(s.window[i])
            else:
                if i < 0:
                    ext.append(self.DL.bit(i + s.phase_L))
                else:
                    ext.append(self.DR.bit((i - W) + s.phase_R))

        # apply rule to get new window
        new = []
        for i in range(1, W+1):
            new.append(self.f(ext[i-1], ext[i], ext[i+1]))
        new = tuple(new)

        # --- keep your heuristic shift, but we will VERIFY later ---
        shift = 0
        center = mid
        if new[center] == self.DL.bit(0):
            shift = -1
        elif new[center] == self.DR.bit(0):
            shift = +1

        new_phase_L = (s.phase_L + shift) % self.DL.P
        new_phase_R = (s.phase_R + shift) % self.DR.P

        return Edge(State(new_phase_L, new_phase_R, new), shift)

    def build(self, max_states: int = 200000):
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

    def find_one_cycle(self, max_len: int = 60) -> Optional[Tuple[List[State], List[int]]]:
        # returns (states along cycle including start repeated at end, shifts along edges)
        for start in self.states:
            stack = [(start, [start], [])]
            while stack:
                cur, path, shifts = stack.pop()
                if len(path) > max_len:
                    continue
                for e in self.graph.get(cur, []):
                    if e.to == start and len(path) >= 2:
                        return (path + [start], shifts + [e.shift])
                    if e.to in path:
                        continue
                    stack.append((e.to, path + [e.to], shifts + [e.shift]))
        return None

# ---------- Build a witness configuration ----------
def build_witness_on_ring(
    DL: Domain,
    DR: Domain,
    W: int,
    state: State,
    L: int = 200,
    cut: int = None,
) -> List[int]:
    if cut is None:
        cut = L // 2
    mid = W // 2

    cfg = [0] * L
    # background
    for i in range(L):
        if i < cut:
            cfg[i] = DL.bit(i + state.phase_L)
        else:
            cfg[i] = DR.bit((i - cut) + state.phase_R)

    # embed state.window centered at cut
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
    DL = Domain((0,1), "01")
    DR = Domain((0,0,1), "001")
    W = 7
    rule_num = 69

    bm = BoundaryMachine(rule_num, DL, DR, W)
    bm.build()
    cyc = bm.find_one_cycle(max_len=80)

    if cyc is None:
        print("No cycle found.")
        raise SystemExit(0)

    states, shifts = cyc
    p = len(shifts)
    v = sum(shifts)
    print(f"Cycle length p={p}, net shift v={v}, speed={v/p:.3f}")

    # witness from first state
    s0 = states[0]
    cfg0 = build_witness_on_ring(DL, DR, W, s0, L=240)
    ok = verify_shift(rule_num, cfg0, p, v)
    print("VERIFY F^p(cfg0) == shift(cfg0, v):", ok)