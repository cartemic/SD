# -*- coding: utf-8 -*-
"""
Functions for finding state properties in frozen and equilibrium flows.

Original functions from Shock and Detonation Toolbox
http://www.galcit.caltech.edu/EDL/public/cantera/html/SD_Toolbox/
"""


def _estimate_cj(
        initial_temperature,
        initial_pressure
):
    # make some guesses to get the working gas closer to CJ state and cut down
    # on the number of iterations required
    guess_gamma = 1.3
    guess_cj_mach = 10

    guess_temperature = 7 * initial_temperature
    guess_pressure = initial_pressure * \
        (
            (1 + guess_gamma) / guess_gamma
        ) * guess_cj_mach

    return {
        'pressure': guess_pressure,
        'temperature': guess_temperature
    }


def get_equilibrium_properties(gas, density, temperature):
    """
    This function calculates the equilibrium pressure and enthalpy given
    temperature and density

    Original function: eq_state in Thermo.py

    Parameters
    ----------
    gas : cantera.composite.Solution
        Working gas object.
    density : float
        Mixture density in kg/m^3.
    temperature : float
        Mixture temperature in K.

    Returns
    -------
    dict
        Dictionary containing pressure and temperature values
    """

    gas.TD = temperature, density
    gas.equilibrate('TV')
    pressure = gas.P
    enthalpy = gas.enthalpy_mass

    return {'pressure': pressure, 'enthalpy': enthalpy}


def get_frozen_properties(gas, density, temperature):
    """
    This function calculates the frozen pressure and enthalpy given temperature
    and density

    Original function: state in Thermo.py

    Parameters
    ----------
    gas : cantera.composite.Solution
        Working gas object.
    density : float
        Mixture density in kg/m^3.
    temperature : float
        Mixture temperature in K.

    Returns
    -------
    dict
        Dictionary containing pressure and temperature values
    """

    gas.TD = temperature, density
    pressure = gas.P
    enthalpy = gas.enthalpy_mass

    return {'pressure': pressure, 'enthalpy': enthalpy}
