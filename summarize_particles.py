from math import gcd
from direct_search_fast import search_direct_fast

def primitive(p: int, v: int):
    g = gcd(p, abs(v))
    return (p // g, v // g)

def run():
    rule_num = 69
    left_pat = (0, 1)        # 01
    right_pat = (0, 0, 1)    # 001

    res = search_direct_fast(
        rule_num, left_pat, right_pat,
        L=240, cut=120,
        max_defect=12,
        max_p=120,
        margin=6,
    )

    # group by primitive (p,v)
    groups = {}
    for (p, v), defect in res.items():
        pp, vv = primitive(p, v)
        groups.setdefault((pp, vv), []).append((p, v, defect))

    # print summary
    print(f"raw pairs: {len(res)}")
    print(f"primitive groups: {len(groups)}\n")

    for (pp, vv) in sorted(groups.keys(), key=lambda x: (x[0], abs(x[1]), x[1])):
        examples = sorted(groups[(pp, vv)], key=lambda t: (t[0], abs(t[1])))
        p, v, d = examples[0]
        dstr = ''.join(map(str, d))
        print(f"primitive p={pp:3d} v={vv:4d} speed={vv/pp: .6f}  example: (p={p},v={v}) defect={dstr}  (#={len(examples)})")

if __name__ == "__main__":
    run()