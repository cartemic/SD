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
from . import tools, states, calculate_error


def get_reflected_frozen_state_0(
        initial_state_gas,
        post_shock_gas,
        working_gas,
        incident_shock_speed):
    """

    reflected_fr
    Calculates frozen post-reflected-shock state assumming u1 = 0

    FUNCTION
    SYNTAX
    [p3,UR,working_gas] = reflected_fr(gas1,post_shock_gas,working_gas,UI)

    INPUT:
    gas1 = gas object at initial state
    post_shock_gas = gas object at post-incident-shock state (already computed)
    working_gas = working gas object
    UI = incident shock speed (m/s)

    OUTPUT:
    p3 = post-reflected-shock pressure (Pa)
    UR = reflected shock speed (m/s)
    working_gas = gas object at frozen post-reflected-shock state

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


def get_reflected_frozen_state(
        particle_speed,
        post_shock_gas,
        working_gas,
        max_iterations=500):
    """

    PostReflectedShock_fr
    Calculates frozen post-Reflected-shock state for a specified shock velocity

    FUNCTION
    SYNTAX
    [working_gas] = PostReflectedShock_fr(u2,post_shock_gas,working_gas)

    INPUT
    u2 = current post-incident-shock lab frame particle speed
    post_shock_gas = gas object at post-incident-shock state (already computed)
    working_gas = working gas object

    OUTPUT
    working_gas = gas object at frozen post-reflected-shock state

    """
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


def get_reflected_equil_state_0(
        initial_state_gas,
        post_shock_gas,
        working_gas,
        incident_shock_speed):
    """

    reflected_eq
    Calculates equilibrium post-reflected-shock state assumming u1 = 0

    FUNCTION
    SYNTAX:
    [p3,UR,gas3] = reflected_eq(gas1,gas2,gas3,UI)

    INPUT:
    gas1 = gas object at initial state
    gas2 = gas object at post-incident-shock state (already computed)
    gas3 = working gas object
    UI = incident shock speed (m/s)

    OUTPUT:
    p3 = post-reflected-shock pressure (Pa)
    UR = reflected shock speed (m/s)
    gas3 = gas object at equilibrium post-reflected-shock state

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
        max_iterations=500):
    """

    PostReflectedShock_eq
    Calculates equilibrium post-reflected-shock state for a specified shock
    velocity

    FUNCTION
    SYNTAX
    [gas3] = PostReflectedShock_fr(u2,gas2,gas3)

    INPUT
    u2 = current post-incident-shock lab frame particle speed
    gas2 = gas object at post-incident-shock state (already computed)
    gas3 = working gas object

    OUTPUT
    gas3 = gas object at equilibrium post-reflected-shock state

    """
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
            states.get_equilibrium_properties(
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
            states.get_equilibrium_properties(
                working_gas,
                guess['density'],
                guess['temperature']
                )
            )

    return working_gas
