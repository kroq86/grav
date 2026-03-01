import unittest

from grav.eca import int_to_list, list_to_int, step_bitwise, step_list


class TestEcaStep(unittest.TestCase):
    def test_bitwise_matches_list_for_all_rules_on_small_ring(self) -> None:
        size = 5
        for rule in range(256):
            for x in range(1 << size):
                state = int_to_list(x, size)
                expected = step_list(state, rule)
                observed = int_to_list(step_bitwise(x, size, rule), size)
                self.assertEqual(expected, observed, msg=f"rule={rule}, state={state}")

    def test_rule_zero_goes_to_all_zero(self) -> None:
        self.assertEqual(step_list([1, 0, 1, 1, 0], 0), [0, 0, 0, 0, 0])

    def test_roundtrip_int_encoding(self) -> None:
        state = [1, 0, 1, 1, 0, 0, 1]
        self.assertEqual(int_to_list(list_to_int(state), len(state)), state)
