# -*- coding: utf-8 -*-
"""
muh docstring
"""

import warnings
import numpy as np
from . import tools, calculate_error


def using_reynolds_iterative_method(working_gas,
                                    initial_state_gas,
                                    error_tol_temperature,
                                    error_tol_specific_volume,
                                    density_ratio,
                                    max_iterations=500):
    """
    This function calculates the Chapman-Jouguet wave speed using Reynolds'
    iterative method.

    Original function: CJ_calc in PostShock.py

    Parameters
    ----------
    working_gas : cantera gas object
        A cantera gas object used for calculations (???).
    initial_state_gas : cantera gas object
        A cantera gas object for the working gas mixture in its initial,
        undetonated state.
    error_tol_temperature : float
        Temperature error tolerance for iteration.
    error_tol_specific_volume : float
        Specific volume error tolerance for iteration.
    density_ratio : float
        density ratio (???).
    max_iterations : int
        Maximum number of loop iterations used to calculate output. Default
        is 500.

    Returns
    -------
    working_gas : cantera gas object
        Gas object at equilibrium state.
    initial_velocity : float
        Initial velocity resulting in the input density ratio, in m/s.
    """
    max_iterations = int(max_iterations)

    # Set initial state
    initial = {
               'density': initial_state_gas.density,
               'volume': 1 / initial_state_gas.density,
               'pressure': initial_state_gas.P
               }

    # Set initial delta values
    delta = {
             'temperature': 1000.,
             'volume': 1000.,
             'pressure': 1000.,
             'velocity': 1000.
             }

    # Set guess values
    guess = {
             'temperature': 2000.,
             'volume': initial['volume'] / density_ratio,
             'density': density_ratio / initial['volume'],
             'initial_velocity': 2000.
             }

    # Add pressure and enthalpy to guess dictionary
    guess.update(eq_state(working_gas,
                          guess['density'],
                          guess['temperature']))

    # NOTE: make eq_state output a dictionary with values for 'pressure' and
    # 'enthalpy'

    # START LOOP
    loop_counter = 0
    while ((abs(delta['temperature']) >
           error_tol_temperature * guess['temperature'])
           or
           (abs(delta['velocity']) >
           error_tol_specific_volume * guess['initial_velocity'])):

        # Manage loop count
        loop_counter += 1
        if loop_counter == max_iterations:
            warnings.warn('No convergence within {0} iterations'.
                          format(max_iterations))
            return [None, None]

        # Calculate pressure and enthalpy error for current guess values as
        # a numpy array, which will be the negative of the ordinate vector, b,
        # for the linear equation Ax = b
        ordinate_vector = (
                           calculate_error.
                           equilibrium(working_gas,
                                       initial_state_gas,
                                       guess['initial_velocity']) * -1
                           )

        delta['temperature'] = guess['temperature'] * 0.02

        # Perturb temperature guess and find corresponding pressure and
        # enthalpy values
        perturbed = tools.perturb('temperature', guess, delta)
        perturbed.update(eq_state(working_gas,
                                  perturbed['density'],
                                  perturbed['temperature']))

        # Calculate pressure and enthalpy error for temperature-perturbed
        # state, add (negative) ordinate vector to get error deltas for
        # enthalpy and pressure, and divide by temperature delta to get rates
        # WRT temperature
        temperature_partials = (
                                calculate_error.
                                equilibrium(
                                            working_gas,
                                            initial_state_gas,
                                            perturbed['temperature']
                                            ) +
                                ordinate_vector
                                )
        temperature_partials /= delta['temperature']

        # Perturb velocity guess -- temperature and enthalpy are not affected
        delta['velocity'] = guess['initial_velocity'] * 0.02
        perturbed = tools.perturb('initial_velocity', guess, delta)

        # Calculate pressure and enthalpy error for velocity-perturbed
        # state, add (negative) ordinate vector to get error deltas for
        # enthalpy and pressure, and divide by velocity delta to get rates
        # WRT velocity
        velocity_partials = (
                             calculate_error.
                             equilibrium(
                                         working_gas,
                                         initial_state_gas,
                                         perturbed['initial_velocity']
                                         ) +
                             ordinate_vector
                             )
        velocity_partials /= delta['initial_velocity']

        # build coefficient matrix, A, for Ax = b
        coefficient_matrix = np.array(temperature_partials,
                                      velocity_partials).transpose()

        # Solve Ax = b
        [delta['temperature'],
         delta['initial_velocity']] = np.linalg.solve(coefficient_matrix,
                                                      ordinate_vector)

        # CHECK & LIMIT CHANGE VALUES
        # VOLUME
        DTM = guess['temperature'] * 0.2
        if abs(delta['temperature']) > DTM:
            delta['temperature'] *= DTM / abs(delta['temperature'])

        # Update guess values
        guess['temperature'] += delta['temperature']
        guess['initial_velocity'] += delta['initial_velocity']
        [P, H] = eq_state(working_gas, guess['density'], guess['temperature'])

    return [working_gas, guess['initial_velocity']]
