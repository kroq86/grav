from __future__ import annotations
from itertools import product
from typing import Tuple, List, Optional

from grav.research.find_translating_fronts import find_translating_fronts

def main():
    rule_num = 69
    domains = [
        ("01", (0, 1)),
        ("001", (0, 0, 1)),
        ("00101", (0, 0, 1, 0, 1)),
        ("0010101", (0, 0, 1, 0, 1, 0, 1)),
        ("00100101", (0, 0, 1, 0, 0, 1, 0, 1)),
    ]

    L = 240
    cut = 120
    max_defect = 12
    margin = 6

    rows = []
    for i in range(len(domains)):
        for j in range(len(domains)):
            if i == j:
                continue
            nameL, patL = domains[i]
            nameR, patR = domains[j]
            res = find_translating_fronts(
                rule_num, patL, patR,
                shifts=(-1, +1),
                L=L, cut=cut,
                max_defect=max_defect,
                margin=margin,
                max_witness_per_shift=1,
            )
            for s in (-1, +1):
                if res[s]:
                    dstr = "".join(map(str, res[s][0]))
                    rows.append((nameL, nameR, s, dstr))

    # print markdown
    print("| Left domain | Right domain | shift per step | witness defect |")
    print("|---|---|---:|---|")
    for a, b, s, d in rows:
        print(f"| `{a}` | `{b}` | {s:+d} | `{d}` |")

if __name__ == "__main__":
    main()
