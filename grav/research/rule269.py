from __future__ import annotations
from dataclasses import dataclass
from itertools import product
from typing import Tuple, Dict, List, Set

# ---------- ECA ----------
def rule_from_number(n: int):
    def f(a: int, b: int, c: int) -> int:
        idx = (1 - a) * 4 + (1 - b) * 2 + (1 - c)  # 111->0 ... 000->7
        return (n >> (7 - idx)) & 1
    return f

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
    shift: int  # -1, 0, +1 per step (radius=1 bound)

class BoundaryMachine:
    def __init__(self, rule_num: int, DL: Domain, DR: Domain, W: int):
        self.f = rule_from_number(rule_num)
        self.DL = DL
        self.DR = DR
        self.W = W  # window size (odd recommended)
        self.states: Set[State] = set()
        self.graph: Dict[State, List[Edge]] = {}

    def initial_states(self) -> List[State]:
        # all possible windows that differ from pure background near interface
        W = self.W
        inits = []
        for w in product([0,1], repeat=W):
            # avoid trivial background windows (both sides match domains exactly)
            # build background window for phase 0|0
            bg = []
            mid = W // 2
            for i in range(W):
                if i < mid:
                    bg.append(self.DL.bit(i - mid))
                else:
                    bg.append(self.DR.bit(i - mid))
            if tuple(bg) == w:
                continue
            inits.append(State(0, 0, w))
        return inits

    def step_state(self, s: State) -> Edge:
        W = self.W
        mid = W // 2

        # reconstruct extended local config of length W+2
        ext = []
        for i in range(-1, W+1):
            j = i
            if 0 <= j < W:
                ext.append(s.window[j])
            else:
                # outside window: use domain bits with current phases
                if j < 0:
                    ext.append(self.DL.bit(j + s.phase_L))
                else:
                    ext.append(self.DR.bit(j - W + s.phase_R))
        # apply rule
        new = []
        for i in range(1, W+1):
            new.append(self.f(ext[i-1], ext[i], ext[i+1]))
        new = tuple(new)

        # estimate shift by comparing center bit with domains
        # if center now matches left/right better, move boundary
        shift = 0
        center = mid
        if new[center] == self.DL.bit(0):
            shift = -1
        elif new[center] == self.DR.bit(0):
            shift = +1

        # update phases
        new_phase_L = (s.phase_L + shift) % self.DL.P
        new_phase_R = (s.phase_R + shift) % self.DR.P

        return Edge(State(new_phase_L, new_phase_R, new), shift)

    def build(self, max_states: int = 50000):
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

    def find_cycles(self, max_len: int = 50):
        cycles = []
        visited = set()

        def dfs(start, cur, path, shifts):
            if len(path) > max_len:
                return
            for e in self.graph.get(cur, []):
                if e.to == start:
                    p = len(path) + 1
                    v = sum(shifts) + e.shift
                    if v != 0:  # moving particles only
                        cycles.append((p, v))
                if e.to not in path:
                    dfs(start, e.to, path | {e.to}, shifts + [e.shift])

        for s in list(self.states):
            dfs(s, s, {s}, [])
        return cycles

# ---------- Run ----------
if __name__ == "__main__":
    DL = Domain((0,1), "01")
    DR = Domain((0,0,1), "001")

    bm = BoundaryMachine(rule_num=69, DL=DL, DR=DR, W=7)
    bm.build()
    cycles = bm.find_cycles(max_len=30)

    # unique (p,v)
    uniq = sorted(set(cycles), key=lambda x: (x[0], abs(x[1]), x[1]))
    print("Found moving cycles (p,v):")
    for p, v in uniq[:20]:
        print(f"  p={p:2d} v={v:3d} speed={v/p: .3f}")