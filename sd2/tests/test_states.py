# -*- coding: utf-8 -*-
"""
Test functions in states.py
"""

import cantera as ct
from .. import states


def test_equilibrium():
    """
    Test the equilibrium state calculator
    """
    # Initialize gas object and collect function result
    test_gas = ct.Solution('gri30.cti')
    test_gas.TPX = 300, 101325, 'H2:2, O2:1'
    test_result = states.get_equilibrium_properties(
        test_gas,
        test_gas.density,
        test_gas.T
        )

    # Equilibrate holding temperature and volume constant, and make sure the
    # answer is as expected
    test_gas.equilibrate('TV')
    [test_pressure, test_enthalpy] = [test_gas.P, test_gas.enthalpy_mass]
    assert abs(test_pressure - test_result['pressure']) <= 1e-7
    assert abs(test_enthalpy - test_result['enthalpy']) <= 1e-7


def test_frozen():
    """
    Test whether or not you want to build a snowman
    """
    # Initialize gas object and collect function result
    test_gas = ct.Solution('gri30.cti')
    test_gas.TPX = 300, 101325, 'H2:2, O2:1'
    test_result = states.get_frozen_properties(
        test_gas,
        test_gas.density,
        test_gas.T
        )

    # Make sure answer matches
    [test_pressure, test_enthalpy] = [test_gas.P, test_gas.enthalpy_mass]
    assert abs(test_pressure - test_result['pressure']) <= 1e-7
    assert abs(test_enthalpy - test_result['enthalpy']) <= 1e-7
