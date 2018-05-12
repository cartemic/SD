# -*- coding: utf-8 -*-
"""
Functions for calculating error in enthalpy and pressure for frozen and
equilibrium gas states.

Original functions from Shock and Detonation Toolbox
http://www.galcit.caltech.edu/EDL/public/cantera/html/SD_Toolbox/
"""

import numpy as np


def equilibrium(working_gas,
                initial_state_gas,
                float initial_velocity):
    """
    This function uses the momentum and energy conservation equations to
    calculate error in current pressure and enthalpy guesses. In this case,
    working state is in equilibrium.

    Original function: FHFP_CJ in PostShock.py

    Parameters
    ----------
    working_gas : cantera gas object
        A cantera gas object used for calculations (???).
    initial_state_gas : cantera gas object
        A cantera gas object for the working gas mixture in its initial,
        undetonated state.
    initial_velocity : float
        A guess for the initial velocity in m/s

    Returns
    -------
    numpy array
        A numpy array of errors in [enthalpy, pressure]
    """
    # define initial state
    cdef float initial_pressure = initial_state_gas.P
    cdef float initial_enthalpy = initial_state_gas.enthalpy_mass
    cdef float initial_density = initial_state_gas.density

    # define equilibrium state
    cdef float equil_pressure = working_gas.P
    cdef float equil_enthalpy = working_gas.enthalpy_mass
    cdef float equil_density = working_gas.density
    cdef float equil_velocity = (
            initial_velocity * (initial_density / equil_density)
            )

    initial_v2 = initial_velocity**2
    equil_v2 = equil_velocity**2
    
    enthalpy_error = ((equil_enthalpy +
                       0.5 * equil_v2) -
                      (initial_enthalpy+
                       0.5 * initial_v2))

    pressure_error = ((equil_pressure +
                       equil_density *
                       equil_v2) -
                      (initial_pressure +
                       initial_density *
                       initial_v2))

    return np.array([enthalpy_error, pressure_error])

