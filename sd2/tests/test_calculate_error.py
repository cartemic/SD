#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test functions in calculate_error.py
"""
from cantera import Solution
from .. import calculate_error


def test_equilibrium():
    mechanism = 'gri30.cti'

    # define states
    initial_gas = Solution(mechanism)
    initial_gas.TPX = 300, 101325, {'H2': 1}

    working_gas = Solution(mechanism)
    working_gas.TPX = 300, 101325*2, {'H2': 1}

    velocity_guess = 7

    # errors from hand calculations
    good_errors = [
        -18.375000000,  # enthalpy
        101322.993718076,  # pressure
    ]

    test_errors = calculate_error.equilibrium(
        working_gas=working_gas,
        initial_state_gas=initial_gas,
        initial_velocity_guess=velocity_guess
    )

    for test, good in zip(test_errors, good_errors):
        assert abs(test - good) / good < 1e-7


def test_frozen():
    mechanism = 'gri30.cti'

    # define states
    initial_gas = Solution(mechanism)
    initial_gas.TPX = 300, 101325, {'H2': 1}

    working_gas = Solution(mechanism)
    working_gas.TPX = 300, 101325*2, {'H2': 1}

    velocity_guess = 7

    # errors from hand calculations
    good_errors = [
        -18.375000000,  # enthalpy
        101322.993718076,  # pressure
    ]

    test_errors = calculate_error.frozen(
        working_gas=working_gas,
        initial_state_gas=initial_gas,
        initial_velocity_guess=velocity_guess
    )

    for test, good in zip(test_errors, good_errors):
        assert abs(test - good) / good < 1e-7


def test_reflected_shock_frozen():
    mechanism = 'gri30.cti'

    # define states
    initial_gas = Solution(mechanism)
    initial_gas.TPX = 300, 101325, {'H2': 1}

    working_gas = Solution(mechanism)
    working_gas.TPX = 300, 101325*2, {'H2': 1}

    velocity_guess = 7

    # errors from hand calculations
    good_errors = [
        -73.5,  # enthalpy
        101316.97487230592  # pressure
    ]

    test_errors = calculate_error.reflected_shock_frozen(
        working_gas=working_gas,
        post_shock_gas=initial_gas,
        shock_speed=velocity_guess
    )

    for test, good in zip(test_errors, good_errors):
        assert abs(test - good) / good < 1e-7
