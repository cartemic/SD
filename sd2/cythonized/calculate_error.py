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
                initial_velocity_guess):
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
    initial_velocity_guess : float
        A guess for the initial velocity in m/s

    Returns
    -------
    numpy array
        A numpy array of errors in [enthalpy, pressure]
    """

    initial = {
               'pressure': initial_state_gas.P,
               'enthalpy': initial_state_gas.enthalpy_mass,
               'density': initial_state_gas.density,
               'velocity': initial_velocity_guess
               }

    equilibrium = {
                   'pressure': working_gas.P,
                   'enthalpy': working_gas.enthalpy_mass,
                   'density': working_gas.density,
                   }

    equilibrium['velocity'] = initial['velocity'] * (initial['density'] /
                                                     equilibrium['density'])

    squared_velocity = {
                    'initial': initial['velocity']**2,
                    'equilibrium': equilibrium['velocity']**2
                    }

    enthalpy_error = ((equilibrium['enthalpy'] +
                       0.5 * squared_velocity['equilibrium']) -
                      (initial['enthalpy'] +
                       0.5 * squared_velocity['initial']))

    pressure_error = ((equilibrium['pressure'] +
                       equilibrium['density'] *
                       squared_velocity['equilibrium']) -
                      (initial['pressure'] +
                       initial['density'] *
                       squared_velocity['initial']))

    return np.array([enthalpy_error, pressure_error])


def frozen(working_gas,
           initial_state_gas,
           initial_velocity_guess):
    """
    This function uses the momentum and energy conservation equations to
    calculate error in current pressure and enthalpy guesses. In this case,
    working state is frozen.

    Original function: FHFP_CJ in PostShock.py

    Do you want to build a snowman?

    Parameters
    ----------
    working_gas : cantera gas object
        A cantera gas object used for calculations.
    initial_state_gas : cantera gas object
        A cantera gas object for the working gas mixture in its initial,
        undetonated state.
    initial_velocity_guess : float
        A guess for the initial velocity in m/s

    Returns
    -------
    numpy array
        A numpy array of errors in [enthalpy, pressure]
    """

    initial = {
               'pressure': initial_state_gas.P,
               'enthalpy': initial_state_gas.enthalpy_mass,
               'density': initial_state_gas.density,
               'velocity': initial_velocity_guess
               }

    working = {
              'pressure': working_gas.P,
              'enthalpy': working_gas.enthalpy_mass,
              'density': working_gas.density,
              }

    working['velocity'] = initial['velocity'] * (initial['density'] /
                                                 working['density'])

    # The working squared velocity is the only part of this which differs
    # from the equilibrium() function
    squared_velocity = {
                    'initial': initial['velocity']**2,
                    'working': (
                                    working['velocity'] *
                                    initial['density'] /
                                    working['density']
                                    )**2
                    }

    enthalpy_error = ((working['enthalpy'] +
                       0.5 * squared_velocity['working']) -
                      (initial['enthalpy'] +
                       0.5 * squared_velocity['initial']))

    pressure_error = ((working['pressure'] +
                       working['density'] *
                       squared_velocity['working']) -
                      (initial['pressure'] +
                       initial['density'] *
                       squared_velocity['initial']))

    return np.array([enthalpy_error, pressure_error])


def reflected_shock_frozen(shock_speed,
                           working_gas,
                           post_shock_gas):
    """
    This function uses the momentum and energy conservation equations to
    calculate error in current pressure and enthalpy guesses during reflected
    shock calculations. In this case, working state is frozen.

    Original function: FHFP_reflected_fr in reflections.py

    Parameters
    ----------
    shock_speed : float
        Current post-incident-shock lab frame particle speed
    working_gas : cantera gas object
        A cantera gas object used for calculations.
    post_shock_gas : cantera gas object
        A cantera gas object at post-incident-shock state (already computed)

    Returns
    -------
    numpy array
        A numpy array of errors in [enthalpy, pressure]
    """
    post_shock = {
        'pressure': post_shock_gas.P,
        'enthalpy': post_shock_gas.enthalpy_mass,
        'density': post_shock_gas.density
        }

    working = {
        'pressure': working_gas.P,
        'enthalpy': working_gas.enthalpy_mass,
        'density': working_gas.density
        }

    enthalpy_error = (
        working['enthalpy'] -
        post_shock['enthalpy'] -
        0.5 * (shock_speed**2)*(
            (working['density'] / post_shock['density']) + 1
            ) /
        (working['density'] / post_shock['density'] - 1)
        )

    pressure_error = (
        working['pressure'] -
        post_shock['pressure'] -
        working['density'] *
        (shock_speed**2) / (
            working['density'] /
            post_shock['density']-1
            )
        )

    return [enthalpy_error, pressure_error]
