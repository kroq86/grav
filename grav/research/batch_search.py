from grav.research.direct_search_fast import search_direct_fast

def run():
    rule_num = 69

    left_pat = (0, 1)        # 01
    right_pat = (0, 0, 1)    # 001

    L_values = [240, 360]
    max_defects = [10, 12]
    max_ps = [80, 120]
    margins = [6, 8]

    for L in L_values:
        cut = L // 2
        for md in max_defects:
            for mp in max_ps:
                for margin in margins:
                    res = search_direct_fast(
                        rule_num,
                        left_pat,
                        right_pat,
                        L=L,
                        cut=cut,
                        max_defect=md,
                        max_p=mp,
                        margin=margin
                    )
                    print(
                        f"L={L:3d} md={md:2d} mp={mp:3d} margin={margin:2d}  "
                        f"->  pairs={len(res)}"
                    )

if __name__ == "__main__":
    run()
