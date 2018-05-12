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
                       float error_tol_temperature,
                       float error_tol_velocity,
                       float density_ratio,
                       int max_iterations=500):
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
    error_tol_velocity : float
        Velocity error tolerance for iteration.
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
    
    # Set initial state
    cdef float initial_density = initial_state_gas.density
    cdef float initial_volume = 1 / initial_density
    cdef float initial_pressure = initial_state_gas.P

    # Set initial delta values
    cdef float delta_temperature = 1000
    cdef float delta_volume = 1000
    cdef float delta_pressure = 1000
    cdef float delta_velocity = 1000

    # Set guess values
    cdef float guess_temperature = 2000
    cdef float guess_density = density_ratio / initial_volume
    cdef float guess_velocity = 2000

    # Add pressure and enthalpy to guess dictionary
    cdef float guess_pressure = 2
    cdef float guess_enthalpy = 2
    [guess_pressure, guess_enthalpy] = states.get_equilibrium_properties(
        working_gas,
        guess_density,
        guess_temperature
        )

    cdef float perturbed_temperature = guess_temperature
    cdef float perturbed_pressure = guess_pressure
    cdef float perturbed_enthalpy = guess_enthalpy
    cdef float perturbed_velocity = guess_velocity
    cdef float max_temperature_delta = 2000
    cdef int loop_counter = 0
    while ((abs(delta_temperature) >
            error_tol_temperature * guess_temperature)
           or
           (abs(delta_velocity) >
            error_tol_velocity * guess_velocity)):

        # Manage loop count
        loop_counter += 1
        if loop_counter == max_iterations:
            warnings.warn('No convergence within {0} iterations'.
                          format(max_iterations))
            return [working_gas, guess_velocity]

        # Calculate pressure and enthalpy error for current guess values as
        # a numpy array, which will be the negative of the ordinate vector, b,
        # for the linear equation Ax = b
        ordinate_vector = (
            calculate_error.
            equilibrium(
                working_gas,
                initial_state_gas,
                guess_velocity) * -1
            )

        # Perturb temperature guess and find corresponding pressure and
        # enthalpy values
        delta_temperature = guess_temperature * 0.02
        perturbed_temperature = guess_temperature + delta_temperature
        [pertubed_pressure,
         perturbed_enthalpy] = states.get_equilibrium_properties(
                working_gas,
                guess_density,
                perturbed_temperature
                )
        # working_gas.TP = perturbed_temperature, perturbed_pressure

        # Calculate pressure and enthalpy error for temperature-perturbed
        # state, add (negative) ordinate vector to get error deltas for
        # enthalpy and pressure, and divide by temperature delta to get rates
        # WRT temperature
        temperature_partials = (
            calculate_error.
            equilibrium(
                working_gas,
                initial_state_gas,
                guess_velocity
                ) +
            ordinate_vector
            )
        temperature_partials /= delta_temperature

        # Perturb velocity guess -- temperature and enthalpy are not affected
        delta_velocity = guess_velocity * 0.02
        perturbed_velocity = guess_velocity + delta_velocity
        # working_gas.TP = guess_temperature, guess_pressure

        # Calculate pressure and enthalpy error for velocity-perturbed
        # state, add (negative) ordinate vector to get error deltas for
        # enthalpy and pressure, and divide by velocity delta to get rates
        # WRT velocity
        velocity_partials = (
            calculate_error.
            equilibrium(
                working_gas,
                initial_state_gas,
                perturbed_velocity
                ) +
            ordinate_vector
            )
        velocity_partials /= delta_velocity

        # build coefficient matrix, A, for Ax = b
        coefficient_matrix = np.array([temperature_partials,
                                       velocity_partials]).transpose()

        # Solve Ax = b
        [delta_temperature,
         delta_velocity] = np.linalg.solve(
             coefficient_matrix,
             ordinate_vector
             )

        # Limit temperature delta to 20% of current guess
        max_temperature_delta = guess_temperature * 0.2
        if abs(delta_temperature) > max_temperature_delta:
            delta_temperature *= (
                max_temperature_delta / abs(delta_temperature)
                )

        # Update guess values
        guess_temperature += delta_temperature
        guess_velocity += delta_velocity
        [guess_pressure,
         guess_enthalpy] = states.get_equilibrium_properties(
                working_gas,
                guess_density,
                guess_temperature
                )

    return [working_gas, guess_velocity]


def calculate_cj_speed(float initial_pressure,
                       float initial_temperature,
                       str species_mole_fractions,
                       str mechanism,
                       return_r_squared=False):
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
    cdef int numsteps = 20
    cdef float max_density_ratio = 2.0
    cdef float min_density_ratio = 1.5

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
    cj_speed = 0.0
    adjusted_density_ratio = 0.0

    def curve_fit_function(x, a, b, c):
        """
        Quadratic function for least-squares curve fit of cj speed vs. density
        ratio
        """
        return a * x**2 + b * x + c

    while (counter <= 4) or (r_squared < 0.99999):
        step = (max_density_ratio - min_density_ratio) / float(numsteps)
        states_calculated = 0
        density_ratio = min_density_ratio
        while density_ratio <= max_density_ratio:
            working_gas.TPX = [
                initial_temperature,
                initial_pressure,
                species_mole_fractions
                ]
            [working_gas,
             temp] = calculate_cj_state(
                 working_gas,
                 initial_state_gas,
                 error_tol_temperature,
                 error_tol_specific_volume,
                 density_ratio
                 )
            cj_velocity_calculations[states_calculated] = temp
            density_ratio_calculations[states_calculated] = (
                working_gas.density /
                initial_state_gas.density
                )
            states_calculated += 1
            density_ratio += step

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
        r_squared = 1 - (
            np.sum(residuals**2) /
            np.sum(
                (
                    cj_velocity_calculations -
                    np.mean(cj_velocity_calculations)
                    )**2
                )
            )

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
    if return_r_squared:
        return [cj_speed, r_squared]
    else:
        return cj_speed

