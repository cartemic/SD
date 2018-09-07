# -*- coding: utf-8 -*-
"""
Functions for detonation calculations.

Original functions from Shock and Detonation Toolbox
http://www.galcit.caltech.edu/EDL/public/cantera/html/SD_Toolbox/
"""

import warnings
import numpy as np
import cantera as ct
from scipy.optimize import curve_fit
from . import tools, calculate_error, states


def calculate_cj_state(working_gas,
                       initial_state_gas,
                       error_tol_temperature,
                       error_tol_specific_volume,
                       density_ratio,
                       max_iterations=500):
    """
    This function calculates the Chapman-Jouguet state and wave speed using
    Reynolds' iterative method.

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
        'velocity': 2000.
        }

    # Add pressure and enthalpy to guess dictionary
    guess.update(states.get_equilibrium_properties(
        working_gas,
        guess['density'],
        guess['temperature'])
    )

    loop_counter = 0
    while ((abs(delta['temperature']) >
            error_tol_temperature * guess['temperature'])
           or
           (abs(delta['velocity']) >
            error_tol_specific_volume * guess['velocity'])):

        # Manage loop count
        loop_counter += 1
        if loop_counter == max_iterations:
            warnings.warn('No convergence within {0} iterations'.
                          format(max_iterations))
            return [working_gas, guess['velocity']]

        # Calculate pressure and enthalpy error for current guess values as
        # a numpy array, which will be the negative of the ordinate vector, b,
        # for the linear equation Ax = b
        ordinate_vector = (
            calculate_error.
            equilibrium(
                working_gas,
                initial_state_gas,
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
            equilibrium(
                working_gas,
                initial_state_gas,
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
            equilibrium(
                working_gas,
                initial_state_gas,
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

    return [working_gas, guess['velocity']]


def calculate_cj_speed(initial_pressure,
                       initial_temperature,
                       species_mole_fractions,
                       mechanism,
                       return_r_squared=False,
                       return_state=False):
    """
    This function calculates CJ detonation velocity

    Original function: CJspeed in PostShock.py

    Parameters
    ----------
    initial_pressure : float
        initial pressure (Pa)
    initial_temperature : float
        initial temperature (K)
    species_mole_fractions : str
        string of reactant species mole fractions
    mechanism : str
        cti file containing mechanism data (e.g. 'gri30.cti')

    Returns
    -------
    cj_speed : float
        CJ detonation speed (m/s)
    r_squared : float
        R-squared value of least-squares parabolic curve fit
    """
    # DECLARATIONS
    numsteps = 20
    max_density_ratio = 2.0
    min_density_ratio = 1.5

    cj_velocity_calculations = np.zeros(numsteps+1, float)
    density_ratio_calculations = np.zeros(numsteps+1, float)

    # Initialize gas objects
    initial_state_gas = ct.Solution(mechanism)
    working_gas = ct.Solution(mechanism)
    initial_state_gas.TPX = [
        initial_temperature,
        initial_pressure,
        species_mole_fractions
        ]
    working_gas.TPX = [
        initial_temperature,
        initial_pressure,
        species_mole_fractions
        ]

    # Set error tolerances for CJ state calculation
    error_tol_temperature = 1e-4
    error_tol_specific_volume = 1e-4

    counter = 1
    r_squared = 0.0
    adjusted_density_ratio = 0.0

    def curve_fit_function(x, a, b, c):
        """
        Quadratic function for least-squares curve fit of cj speed vs. density
        ratio
        """
        return a * x**2 + b * x + c

    def calculate_over_ratio_range(
            current_state_number,
            current_density_ratio,
            current_gas
    ):
        current_gas.TPX = [
            initial_temperature,
            initial_pressure,
            species_mole_fractions
        ]
        [current_gas,
         temp] = calculate_cj_state(
            working_gas,
            initial_state_gas,
            error_tol_temperature,
            error_tol_specific_volume,
            current_density_ratio
        )
        current_velocity = temp
        density_ratio_calculations[states_calculated] = (
                current_gas.density /
                initial_state_gas.density
        ) # find out if this is any different brb.


    while (counter <= 4) and (r_squared < 0.99999 or delta_r_squared < 1e-7):
        density_ratio_array = np.linspace(
            min_density_ratio,
            max_density_ratio,
            numsteps + 1
        )
        for states_calculated, density_ratio in enumerate(density_ratio_array):
            # add thing here

        # Get curve fit
        [curve_fit_coefficients, _] = curve_fit(
            curve_fit_function,
            density_ratio_calculations,
            cj_velocity_calculations
            )

        # Calculate R^2 value
        residuals = cj_velocity_calculations - curve_fit_function(
            density_ratio_calculations,
            *curve_fit_coefficients
            )
        old_r_squared = r_squared
        r_squared = 1 - (
            np.sum(residuals**2) /
            np.sum(
                (
                    cj_velocity_calculations -
                    np.mean(cj_velocity_calculations)
                    )**2
                )
            )
        delta_r_squared = abs(old_r_squared - r_squared)

        adjusted_density_ratio = (
            -curve_fit_coefficients[1] /
            (2. * curve_fit_coefficients[0])
            )
        min_density_ratio = adjusted_density_ratio * (1 - 0.001)
        max_density_ratio = adjusted_density_ratio * (1 + 0.001)
        counter += 1

    cj_speed = curve_fit_function(
        adjusted_density_ratio,
        *curve_fit_coefficients
        )

    if return_r_squared and return_state:
        return [cj_speed, r_squared, working_gas]
    elif return_state:
        return [cj_speed, working_gas]
    elif return_r_squared:
        return [cj_speed, r_squared]
    else:
        return cj_speed
