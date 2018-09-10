#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test functions in detonations.py
"""

from .. import detonations


def test_calculate_cj_speed():
    """
    Test calculate_cj.speed() against known values from original SDToolbox
    """
    test_args = [
        101325,
        300,
        'H2:0.5333 O2:0.26667 AR:0.2', 'gri30.cti'
        ]

    # CJ speed calculated by SDToolbox given *test_args
    original_cj_speed = 2353.2706464533471

    test_speed = detonations.calculate_cj_speed(*test_args)['cj speed']
    relative_difference = abs(
        (original_cj_speed - test_speed) /
        original_cj_speed
        )
    assert relative_difference <= 1e-5
