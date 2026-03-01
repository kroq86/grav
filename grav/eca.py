from __future__ import annotations


def rule_from_number(rule_num: int):
    """Return a Wolfram ECA local update function for triples (left, center, right)."""
    if not 0 <= rule_num <= 255:
        raise ValueError("rule_num must be in [0, 255]")

    def step_cell(left: int, center: int, right: int) -> int:
        idx = (left << 2) | (center << 1) | right
        return (rule_num >> idx) & 1

    return step_cell


def step_list(state: list[int], rule_num: int) -> list[int]:
    """Advance one time step on a periodic ring."""
    step_cell = rule_from_number(rule_num)
    size = len(state)
    return [
        step_cell(state[(i - 1) % size], state[i], state[(i + 1) % size])
        for i in range(size)
    ]


def mask_L(size: int) -> int:
    return (1 << size) - 1


def rot_r(x: int, size: int, shift: int) -> int:
    shift %= size
    masked = x & mask_L(size)
    if shift == 0:
        return masked
    return ((masked >> shift) | ((masked << (size - shift)) & mask_L(size))) & mask_L(size)


def step_bitwise(x: int, size: int, rule_num: int) -> int:
    """Bitwise ECA step with bit i corresponding to cell i."""
    masked = x & mask_L(size)
    left = ((masked << 1) & mask_L(size)) | (masked >> (size - 1))
    right = (masked >> 1) | ((masked & 1) << (size - 1))
    center = masked

    nxt = 0
    for a in (0, 1):
        left_mask = left if a else ~left
        for b in (0, 1):
            center_mask = center if b else ~center
            for c in (0, 1):
                right_mask = right if c else ~right
                pattern_mask = left_mask & center_mask & right_mask & mask_L(size)
                idx = (a << 2) | (b << 1) | c
                if (rule_num >> idx) & 1:
                    nxt |= pattern_mask
    return nxt & mask_L(size)


def list_to_int(state: list[int]) -> int:
    out = 0
    for i, bit in enumerate(state):
        if bit:
            out |= 1 << i
    return out


def int_to_list(x: int, size: int) -> list[int]:
    return [(x >> i) & 1 for i in range(size)]
