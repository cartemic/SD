#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Functions for shock calculations.

Original functions from Shock and Detonation Toolbox
http://www.galcit.caltech.edu/EDL/public/cantera/html/SD_Toolbox/
"""

from math import sqrt
import warnings
import numpy as np
import cantera as ct
from . import tools, states, calculate_error

''' commented out until dependency is fixed
def get_reflected_frozen_state_0(
        initial_state_gas,
        post_shock_gas,
        working_gas,
        incident_shock_speed):
    """
    This function calculates frozen post-reflected-shock state assuming u1 = 0

    Original function: reflected_fr in reflections.py

    Parameters
    ----------
    initial_state_gas : cantera gas object
        gas object at initial state
    post_shock_gas : cantera gas object
        gas object at post-incident-shock state (already computed)
    working_gas : cantera gas object
        working gas object
    incident_shock_speed : float
        incident shock speed (m/s)

    Returns
    -------
    working['pressure'] : float
        post-reflected-shock pressure (Pa)
    reflected_shock_speed : float
        reflected shock speed (m/s)
    working_gas : cantera gas object
        gas object at frozen post-reflected-shock state
    """
    initial = {
        'pressure': initial_state_gas.P,
        'density': initial_state_gas.density,
        'volume': 1 / initial_state_gas.density
        }

    post_shock = {
        'pressure': post_shock_gas.P,
        'density': post_shock_gas.density,
        'volume': 1 / post_shock_gas.density,
        'temperature': post_shock_gas.T
        }

    particle_speed = sqrt(
        (post_shock['pressure'] - initial['pressure']) *
        (initial['volume'] - post_shock['volume'])
        )

    # BASIC PRELIMINARY GUESS
    working = {
        'volume': 0.2 / post_shock['density']
        }

    working.update({
        'pressure':
        post_shock['pressure'] +
        post_shock['density'] *
        (incident_shock_speed**2) *
        (1 - working['volume'] / post_shock['volume'])
        })

    working.update({
        'temperature':
        post_shock['temperature'] *
        working['pressure'] *
        working['volume'] /
        (post_shock['pressure'] * post_shock['volume'])
        })

    # Update working gas properties
    working_gas.TPX = [
        working['temperature'],
        working['pressure'],
        post_shock_gas.X
        ]
    working_gas = get_reflected_frozen_state(
        particle_speed,
        post_shock_gas,
        working_gas
        )

    working['pressure'] = working_gas.P
    reflected_shock_speed = (
        (working['pressure'] - post_shock['pressure']) /
        particle_speed /
        post_shock['density'] -
        particle_speed
        )

    return [working['pressure'], reflected_shock_speed, working_gas]
'''

# TODO: redo this function without mistaking spec volume for velocity >:(
'''
def get_reflected_frozen_state(
        particle_speed,
        post_shock_gas,
        working_gas,
        max_iterations=500):
    """
    This function calculates frozen post-Reflected-shock state for a specified
    shock velocity

    Original function: PostReflectedShock_fr in reflections.py

    Parameters
    ----------
    particle_speed : float
        current post-incident-shock lab frame particle speed
    post_shock_gas : cantera gas object
        gas object at post-incident-shock state (already computed)
    working_gas : cantera gas object
        working gas object
    max_iterations : int
        maximum number of loop iterations

    Returns
    -------
    working_gas : cantera gas object
        gas object at frozen post-reflected-shock state
    """
    # TODO: redo this function without mistaking spec volume for velocity >:(
    # Set error tolerances for state calculation
    error_tol_temperature = 1e-4
    error_tol_specific_volume = 1e-4

    # Set reflected values
    reflected = {
        'volume': 1 / post_shock_gas.density
        }

    # Set deltas
    delta = {
        'temperature': 1000,
        'volume': 1000
        }

    # Set guess values
    guess = {
        'pressure': working_gas.P,
        'enthalpy': working_gas.enthalpy_mass,
        'temperature': working_gas.T,
        'density': working_gas.density,
        'volume': 1 / working_gas.density
        }

    loop_counter = 0
    while ((abs(delta['temperature']) >
            error_tol_temperature*guess['temperature'])
           or
           (abs(delta['volume']) >
            error_tol_specific_volume*guess['volume'])):

        # Manage loop count
        loop_counter += 1
        if loop_counter == max_iterations:
            warnings.warn('Calculation did not converge for U = {0}'.
                          format(particle_speed))
            return working_gas

        # Calculate pressure and enthalpy error for current guess values as
        # a numpy array, which will be the negative of the ordinate vector, b,
        # for the linear equation Ax = b
        ordinate_vector = (
            calculate_error.
            reflected_shock_frozen(
                working_gas,
                post_shock_gas,
                guess['velocity']) * -1
            )

        # Perturb temperature guess and find corresponding pressure and
        # enthalpy values
        delta['temperature'] = guess['temperature'] * 0.02
        perturbed = tools.perturb('temperature', guess, delta)
        perturbed.update(
            states.get_frozen_properties(
                working_gas,
                perturbed['density'],
                perturbed['temperature']
                )
            )

        # Calculate pressure and enthalpy error for temperature-perturbed
        # state, add (negative) ordinate vector to get error deltas for
        # enthalpy and pressure, and divide by temperature delta to get rates
        # WRT temperature
        temperature_partials = (
            calculate_error.
            reflected_shock_frozen(
                working_gas,
                post_shock_gas,
                perturbed['velocity']
                ) +
            ordinate_vector
            )
        temperature_partials /= delta['temperature']

        # Perturb velocity guess -- temperature and enthalpy are not affected
        delta['velocity'] = guess['velocity'] * 0.02
        perturbed = tools.perturb('velocity', guess, delta)

        # Calculate pressure and enthalpy error for velocity-perturbed
        # state, add (negative) ordinate vector to get error deltas for
        # enthalpy and pressure, and divide by velocity delta to get rates
        # WRT velocity
        velocity_partials = (
            calculate_error.
            reflected_shock_frozen(
                working_gas,
                post_shock_gas,
                perturbed['velocity']
                ) +
            ordinate_vector
            )
        velocity_partials /= delta['velocity']

        # build coefficient matrix, A, for Ax = b
        coefficient_matrix = np.array([temperature_partials,
                                       velocity_partials]).transpose()

        # Solve Ax = b
        [delta['temperature'],
         delta['velocity']] = np.linalg.solve(
             coefficient_matrix,
             ordinate_vector
             )

        # Limit temperature delta to 20% of current guess
        max_temperature_delta = guess['temperature'] * 0.2
        if abs(delta['temperature']) > max_temperature_delta:
            delta['temperature'] *= (
                max_temperature_delta / abs(delta['temperature'])
                )

        # Limit volume delta
        max_volume = guess['volume'] + delta['volume']
        if max_volume > reflected['volume']:
            max_volume_delta = 0.5 * (reflected['volume'] - guess['volume'])
        else:
            max_volume_delta = 0.2 * guess['volume']

        if abs(delta['volume']) > max_volume_delta:
            delta['volume'] = max_volume_delta * np.sign(delta['volume'])

        # Update guess values
        guess['temperature'] += delta['temperature']
        guess['velocity'] += delta['velocity']
        guess.update(
            states.get_frozen_properties(
                working_gas,
                guess['density'],
                guess['temperature']
                )
            )

    return working_gas
'''


def get_reflected_equil_state_0(
        initial_state_gas,
        post_shock_gas,
        working_gas,
        incident_shock_speed):
    """
    This function calculates equilibrium post-reflected-shock state assuming
    u1 = 0

    reflected_eq

    Parameters
    ----------
    initial_state_gas : cantera gas object
        gas object at initial state
    post_shock_gas : cantera gas object
        gas object at post-incident-shock state (already computed)
    working_gas : cantera gas object
        working gas object
    incident_shock_speed : float
        incident shock speed (m/s)

    Returns
    -------
    working['pressure'] : float
        post-reflected-shock pressure (Pa)
    reflected_shock_speed : float
        reflected shock speed (m/s)
    working_gas : cantera gas object
        gas object at equilibrium post-reflected-shock state
    """
    initial = {
        'pressure': initial_state_gas.P,
        'volume': 1 / initial_state_gas.density
        }

    reflected = {
        'pressure': post_shock_gas.P,
        'density': post_shock_gas.density,
        'volume': 1 / post_shock_gas.density,
        'temperature': post_shock_gas.T
        }
    reflected['velocity'] = sqrt(
        (reflected['pressure'] - initial['pressure']) *
        (initial['volume'] - reflected['volume'])
        )

    working = {
        'volume': 0.2 / reflected['density']
        }
    working['pressure'] = (
        reflected['pressure'] +
        reflected['density'] *
        (incident_shock_speed**2) *
        (1 - working['volume'] / reflected['volume'])
        )
    working['temperature'] = (
        reflected['temperature'] *
        working['pressure'] *
        working['volume'] /
        (reflected['pressure'] * reflected['volume'])
        )

    working_gas.TPX = [
        working['temperature'],
        working['pressure'],
        post_shock_gas.X
        ]
    working_gas = get_reflected_equil_state(
        reflected['velocity'],
        post_shock_gas,
        working_gas
        )
    working['pressure'] = working_gas.P
    reflected_shock_speed = (
        (working['pressure'] - reflected['pressure']) /
        reflected['velocity'] /
        reflected['density'] -
        reflected['velocity']
        )

    return [working['pressure'], reflected_shock_speed, working_gas]


def get_reflected_equil_state(
        particle_speed,
        post_shock_gas,
        working_gas,
        error_tol_temperature=1e-4,
        error_tol_specific_volume=1e-4,
        max_iterations=500):
    """
    This function calculates equilibrium post-reflected-shock state for a
    specified shock velocity

    Original function: PostReflectedShock_eq in reflections.py

    Parameters
    ----------
    particle_speed : float
        current post-incident-shock lab frame particle speed
    post_shock_gas : cantera gas object
        gas object at post-incident-shock state (already computed)
    working_gas : cantera gas object
        working gas object
    error_tol_temperature : float
        Temperature error tolerance for iteration.
    error_tol_specific_volume : float
        Specific volume error tolerance for iteration.
    max_iterations : int
        maximum number of loop iterations

    Returns
    -------
    working_gas : cantera gas object
        gas object at equilibrium post-reflected-shock state
    """

    # set post-shocked state
    post_shock = dict()
    post_shock['volume'] = 1 / post_shock_gas.density

    # set reflected guess state
    guess = dict()
    guess['temperature'] = working_gas.T
    guess['density'] = working_gas.density
    guess['volume'] = 1 / guess['density']

    # set initial delta guesses
    delta = dict()
    delta['temperature'] = 1000
    delta['volume'] = 1000

    # equilibrate at guess state
    states.get_equilibrium_properties(
        working_gas,
        guess['density'],
        guess['temperature']
    )

    # calculate reflected state
    loop_counter = 0
    while (
            (
                    abs(delta['temperature'])
                    >
                    error_tol_temperature * guess['temperature'])
            or
            (
                    abs(delta['volume'])
                    >
                    error_tol_specific_volume * guess['volume']
            )
    ):
        loop_counter += 1
        if loop_counter == max_iterations:
            warnings.warn(
                'Calculation did not converge for U = {0:.2f} ' +
                'after {1} iterations'.format(particle_speed, loop_counter))
            return working_gas

        # calculate enthalpy and pressure error for current guess
        # TODO: look into eq function instead of frozen
        [err_enthalpy, err_pressure] = calculate_error.reflected_shock_frozen(
            particle_speed,
            working_gas,
            post_shock_gas
        )

        # equilibrate working gas with perturbed temperature
        delta['temperature'] = guess['temperature'] * 0.02
        # equilibrate temperature perturbed state
        states.get_equilibrium_properties(
            working_gas,
            guess['density'],
            guess['temperature'] + delta['temperature']
        )

        # calculate enthalpy and pressure error for perturbed temperature
        [err_enthalpy_perturbed,
         err_pressure_perturbed] = calculate_error.reflected_shock_frozen(
            particle_speed,
            working_gas,
            post_shock_gas
        )

        # calculate temperature derivatives
        deriv_enthalpy_temperature = (err_enthalpy_perturbed -
                                      err_enthalpy) / delta['temperature']
        deriv_pressure_temperature = (err_pressure_perturbed -
                                      err_pressure) / delta['temperature']

        # equilibrate working gas with perturbed volume
        delta['volume'] = 0.02 * guess['volume']
        # equilibrate volume perturbed state
        states.get_equilibrium_properties(
            working_gas,
            1 / (guess['volume'] + delta['volume']),
            guess['temperature']
        )

        # calculate enthalpy and pressure error for perturbed specific volume
        [err_enthalpy_perturbed,
         err_pressure_perturbed] = calculate_error.reflected_shock_frozen(
            particle_speed,
            working_gas,
            post_shock_gas
        )

        # calculate specific volume derivatives
        deriv_enthalpy_volume = (err_enthalpy_perturbed -
                                 err_enthalpy) / delta['volume']
        deriv_pressure_volume = (err_pressure_perturbed -
                                 err_pressure) / delta['volume']

        # solve matrix for temperature and volume deltas
        jacobian = (
                deriv_enthalpy_temperature *
                deriv_pressure_volume -
                deriv_pressure_temperature *
                deriv_enthalpy_volume
        )
        bb = [
            deriv_pressure_volume,
            -deriv_enthalpy_volume,
            -deriv_pressure_temperature,
            deriv_enthalpy_temperature
        ]
        aa = [-err_enthalpy, -err_pressure]
        delta['temperature'] = (bb[0] * aa[0] +
                                bb[1] * aa[1]) / jacobian
        delta['volume'] = (bb[2] * aa[0] +
                           bb[3] * aa[1]) / jacobian

        # check and limit temperature delta
        delta['temp_max'] = 0.2 * guess['temperature']
        if abs(delta['temperature']) > delta['temp_max']:
            delta['temperature'] = (
                    delta['temp_max'] *
                    delta['temperature'] /
                    abs(delta['temperature'])
            )

        # check and limit specific volume delta
        perturbed_volume = guess['volume'] + delta['volume']
        if perturbed_volume > post_shock['volume']:
            delta['volume_max'] = 0.5 * (post_shock['volume'] - guess['volume'])
        else:
            delta['volume_max'] = 0.2 * guess['volume']

        if abs(delta['volume']) > delta['volume_max']:
            delta['volume'] = (
                    delta['volume_max'] *
                    delta['volume'] /
                    abs(delta['volume'])
            )

        # apply calculated and limited deltas to temperature and spec. volume
        guess['temperature'] += + delta['temperature']
        guess['volume'] += delta['volume']
        guess['density'] = 1 / guess['volume']

        # equilibrate working gas with updated state
        states.get_equilibrium_properties(
            working_gas,
            guess['density'],
            guess['temperature']
        )

    return working_gas


def get_post_shock_eq_state(
        wave_speed,
        initial_pressure,
        initial_temperature,
        reactant_mixture,
        mechanism,
        error_tol_temperature=1e-4,
        error_tol_specific_volume=1e-4,
        max_iterations=500):
    """
    This function calculates equilibrium post-shock state using Reynolds'
    iterative method

    Original functions: shk_eq_calc and PostShock_eq in reflections.py

    Parameters
    ----------
    wave_speed : float
        speed at which the shock is traveling
    initial_pressure : float
        Pressure of initial state mixture (Pa)
    initial_temperature : float
        Temperature of initial state mixture (K)
    reactant_mixture : str or dict
        String or dict of reactant species moles or mole fractions
    mechanism : str
        Mechanism file to use (e.g. 'gri30.cti')
    error_tol_temperature : float
        Temperature error tolerance for iteration.
    error_tol_specific_volume : float
        Specific volume error tolerance for iteration.
    max_iterations : int
        maximum number of loop iterations

    Returns
    -------
    working_gas : cantera gas object
        gas object at equilibrium post-reflected-shock state
    """
    # set gas objects
    initial_state_gas = ct.Solution(mechanism)
    initial_state_gas.TPX = [
        initial_temperature,
        initial_pressure,
        reactant_mixture
    ]
    working_gas = ct.Solution(mechanism)
    working_gas.TPX = [
        initial_temperature,
        initial_pressure,
        reactant_mixture
    ]

    # set initial state variables
    initial = dict()
    initial['density'] = initial_state_gas.density
    initial['volume'] = 1/initial['density']
    initial['pressure'] = initial_pressure
    initial['temperature'] = initial_temperature

    # set initial delta guess
    delta = dict()
    delta['temperature'] = 1000
    delta['volume'] = 1000

    # set guess state variables
    guess = dict()
    guess['volume'] = 0.2 * initial['volume']
    guess['density'] = 1 / guess['volume']
    guess['pressure'] = (
            initial['pressure'] +
            initial['density'] *
            (wave_speed**2) *
            (1 - guess['volume'] / initial['volume'])
    )
    guess['temperature'] = (
            initial['temperature'] *
            guess['pressure'] *
            guess['volume'] /
            (initial['pressure'] * initial['volume'])
    )

    # equilibrate working gas
    states.get_equilibrium_properties(
        working_gas,
        guess['density'],
        guess['temperature']
    )

    # calculate equilibrium state
    loop_counter = 0
    while (
            (
                    abs(delta['temperature'])
                    >
                    error_tol_temperature * guess['temperature']
            )
            or
            (
                    abs(delta['volume'])
                    >
                    error_tol_specific_volume*guess['volume']
            )
    ):
        loop_counter += 1
        if loop_counter == max_iterations:
            warnings.warn(
                "No convergence in {0} iterations".format(loop_counter)
            )
            return working_gas

        # calculate enthalpy and pressure error for current guess
        [err_enthalpy, err_pressure] = calculate_error.equilibrium(
            working_gas,
            initial_state_gas,
            wave_speed
        )

        # equilibrate working gas with perturbed temperature
        delta['temperature'] = 0.02 * guess['temperature']
        states.get_equilibrium_properties(
            working_gas,
            guess['density'],
            guess['temperature'] + delta['temperature']
        )

        # calculate enthalpy and pressure error for perturbed temperature
        [err_enthalpy_perturbed,
         err_pressure_perturbed] = calculate_error.equilibrium(
            working_gas,
            initial_state_gas,
            wave_speed
        )

        # calculate temperature derivatives
        deriv_enthalpy_temperature = (err_enthalpy_perturbed -
                                      err_enthalpy) / delta['temperature']
        deriv_pressure_temperature = (err_pressure_perturbed -
                                      err_pressure) / delta['temperature']

        # equilibrate working gas with perturbed volume
        delta['volume'] = 0.02 * guess['volume']
        states.get_equilibrium_properties(
            working_gas,
            1 / (guess['volume'] + delta['volume']),
            guess['temperature']
        )

        # calculate enthalpy and pressure error for perturbed specific volume
        [err_enthalpy_perturbed,
         err_pressure_perturbed] = calculate_error.equilibrium(
            working_gas,
            initial_state_gas,
            wave_speed
        )

        # calculate specific volume derivatives
        deriv_enthalpy_volume = (err_enthalpy_perturbed -
                                 err_enthalpy) / delta['volume']
        deriv_pressure_volume = (err_pressure_perturbed -
                                 err_pressure) / delta['volume']

        # solve matrix for temperature and volume deltas
        jacobian = (
                deriv_enthalpy_temperature *
                deriv_pressure_volume -
                deriv_pressure_temperature *
                deriv_enthalpy_volume
        )
        bb = [
            deriv_pressure_volume,
            -deriv_enthalpy_volume,
            -deriv_pressure_temperature,
            deriv_enthalpy_temperature
        ]
        aa = [-err_enthalpy, -err_pressure]
        delta['temperature'] = (bb[0] * aa[0] +
                                bb[1] * aa[1]) / jacobian
        delta['volume'] = (bb[2] * aa[0] +
                           bb[3] * aa[1]) / jacobian

        # check and limit temperature delta
        delta['temp_max'] = 0.2 * guess['temperature']
        if abs(delta['temperature']) > delta['temp_max']:
            delta['temperature'] = delta['temp_max'] * \
                                   delta['temperature'] / \
                                   abs(delta['temperature'])

        # check and limit specific volume delta
        perturbed_volume = guess['volume'] + delta['volume']
        if perturbed_volume > initial['volume']:
            delta['volume_max'] = 0.5 * (initial['volume'] - guess['volume'])
        else:
            delta['volume_max'] = 0.2 * guess['volume']
        if abs(delta['volume']) > delta['volume_max']:
            delta['volume'] = delta['volume_max'] * \
                              delta['volume'] / \
                              abs(delta['volume'])

        # apply calculated and limited deltas to temperature and spec. volume
        guess['temperature'] += delta['temperature']
        guess['volume'] += delta['volume']
        guess['density'] = 1/guess['volume']

        # equilibrate working gas with updated state
        states.get_equilibrium_properties(
            working_gas,
            guess['density'],
            guess['temperature']
        )

    return working_gas
