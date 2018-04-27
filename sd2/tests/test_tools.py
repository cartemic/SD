# -*- coding: utf-8 -*-
"""
Test functions in tools.py
"""

from .. import tools
from random import uniform, choice
import string


def test_perturb():
    """
    Test the perturb function
    """
    properties = {}
    deltas = {}
    solutions = {}
    for _ in range(3):
        key = ''.join(choice(string.ascii_letters) for _ in range(16))
        value = uniform(-1e6, 1e6)
        delta = uniform(-1e3, 1e3)

        properties[key] = value
        deltas[key] = delta
        solutions[key] = value + delta

    # Ensure correct output
    for property_string in deltas:
        test_result = tools.perturb(property_string, properties, deltas)
        assert test_result[property_string] == solutions[property_string]
