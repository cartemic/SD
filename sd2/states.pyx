# -*- coding: utf-8 -*-
"""
Functions for finding state properties in frozen and equilibrium flows.

Original functions from Shock and Detonation Toolbox
http://www.galcit.caltech.edu/EDL/public/cantera/html/SD_Toolbox/
"""


def get_equilibrium_properties(gas,
                               float density,
                               float temperature):
    """
    This function calculates the equilibrium pressure and enthalpy given
    temperature and density

    Original function: eq_state in Thermo.py

    Parameters
    ----------
    gas : cantera gas object
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
    cdef float pressure = gas.P
    cdef float enthalpy = gas.enthalpy_mass

    return [pressure, enthalpy]
