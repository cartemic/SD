#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test functions in detonations.py
"""
from cantera import Solution
from .. import detonations


class TestCalculateCJSpeed:
    @staticmethod
    def test_no_multiprocessing():
        # Check for correct calculation when multiprocessing is not used
        test_args = [
            101325,
            300,
            'H2:0.5333 O2:0.26667 AR:0.2',
            'gri30.cti',
            False,  # use multiprocessing
            False,  # return r squared
            False   # return state
            ]

        # CJ speed calculated by SDToolbox given *test_args
        original_cj_speed = 2353.2706464533471

        test_speed = detonations.calculate_cj_speed(*test_args)['cj speed']
        relative_difference = abs(original_cj_speed - test_speed) / \
            original_cj_speed

        assert relative_difference <= 1e-5

    @staticmethod
    def test_multiprocessing():
        # Check for correct calculation when multiprocessing is used
        test_args = [
            101325,
            300,
            'H2:0.5333 O2:0.26667 AR:0.2',
            'gri30.cti',
            True,   # use multiprocessing
            False,  # return r squared
            False   # return state
        ]

        # CJ speed calculated by SDToolbox given *test_args
        original_cj_speed = 2353.2706464533471

        test_speed = detonations.calculate_cj_speed(*test_args)['cj speed']
        relative_difference = abs(original_cj_speed - test_speed) / \
            original_cj_speed

        assert relative_difference <= 1e-5

    @staticmethod
    def test_return_state():
        # Check for correct calculation when multiprocessing is used
        test_args = [
            101325,
            300,
            'H2:0.5333 O2:0.26667 AR:0.2',
            'gri30.cti',
            True,   # use multiprocessing
            False,  # return r squared
            True    # return state
        ]
        test_info = detonations.calculate_cj_speed(*test_args)
        required_keys = {
            'cj speed',
            'cj state'
        }

        check_keys = len(required_keys.difference(test_info.keys())) == 0
        check_state = isinstance(test_info['cj state'], Solution)

        assert all([check_keys, check_state])

    @staticmethod
    def test_return_r_squared():
        # Check for correct calculation when multiprocessing is used
        test_args = [
            101325,
            300,
            'H2:0.5333 O2:0.26667 AR:0.2',
            'gri30.cti',
            True,   # use multiprocessing
            True,   # return r squared
            False   # return state
        ]
        test_info = detonations.calculate_cj_speed(*test_args)
        required_keys = {
            'cj speed',
            'R^2'
        }

        check_keys = len(required_keys.difference(test_info.keys())) == 0
        check_r_squared = test_info['R^2'] >= 0.99

        assert all([check_keys, check_r_squared])

    @staticmethod
    def test_return_state_and_r_squared():
        # Check for correct calculation when multiprocessing is used
        test_args = [
            101325,
            300,
            'H2:0.5333 O2:0.26667 AR:0.2',
            'gri30.cti',
            True,   # use multiprocessing
            True,   # return r squared
            True    # return state
        ]
        test_info = detonations.calculate_cj_speed(*test_args)
        required_keys = {
            'cj speed',
            'cj state',
            'R^2'
        }

        check_keys = len(required_keys.difference(test_info.keys())) == 0
        check_state = isinstance(test_info['cj state'], Solution)
        check_r_squared = test_info['R^2'] >= 0.99

        assert all([check_keys, check_state, check_r_squared])


# TODO: tests for calculate_cj_state (including convergence failure!)
