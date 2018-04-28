# -*- coding: utf-8 -*-
"""
Functions for calculating Chapman-Jouguet wave speeds using various solution
methods.

Original functions from Shock and Detonation Toolbox
http://www.galcit.caltech.edu/EDL/public/cantera/html/SD_Toolbox/
"""

import warnings
import numpy as np
import cantera as ct
from scipy.optimize import curve_fit
from . import tools, calculate_error, get_state


def state_and_speed(working_gas,
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
        'initial_velocity': 1000.
        }

    # Set guess values
    guess = {
        'temperature': 2000.,
        'volume': initial['volume'] / density_ratio,
        'density': density_ratio / initial['volume'],
        'initial_velocity': 2000.
        }

    # Add pressure and enthalpy to guess dictionary
    guess.update(get_state.equilibrium(working_gas,
                                       guess['density'],
                                       guess['temperature']))

    loop_counter = 0
    while ((abs(delta['temperature']) >
            error_tol_temperature * guess['temperature'])
           or
           (abs(delta['initial_velocity']) >
            error_tol_specific_volume * guess['initial_velocity'])):

        # Manage loop count
        loop_counter += 1
        if loop_counter == max_iterations:
            warnings.warn('No convergence within {0} iterations'.
                          format(max_iterations))
            return [working_gas, guess['initial_velocity']]

        # Calculate pressure and enthalpy error for current guess values as
        # a numpy array, which will be the negative of the ordinate vector, b,
        # for the linear equation Ax = b
        ordinate_vector = (
            calculate_error.
            equilibrium(
                working_gas,
                initial_state_gas,
                guess['initial_velocity']) * -1
            )

        delta['temperature'] = guess['temperature'] * 0.02

        # Perturb temperature guess and find corresponding pressure and
        # enthalpy values
        perturbed = tools.perturb('temperature', guess, delta)
        perturbed.update(get_state.equilibrium(working_gas,
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
                perturbed['initial_velocity']
                ) +
            ordinate_vector
            )
        temperature_partials /= delta['temperature']

        # Perturb velocity guess -- temperature and enthalpy are not affected
        delta['initial_velocity'] = guess['initial_velocity'] * 0.02
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
        coefficient_matrix = np.array([temperature_partials,
                                       velocity_partials]).transpose()

        # Solve Ax = b
        [delta['temperature'],
         delta['initial_velocity']] = np.linalg.solve(
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
        guess['initial_velocity'] += delta['initial_velocity']
        guess.update(
            get_state.equilibrium(
                working_gas,
                guess['density'],
                guess['temperature']
                )
            )

    return [working_gas, guess['initial_velocity']]


def speed(initial_pressure,
          initial_temperature,
          species_mole_fractions,
          mechanism):
    """

    CJspeed
    Calculates CJ detonation velocity

    FUNCTION
    SYNTAX
    [cj_speed,R2] = CJspeed(P1,T1,q,mech,plt_num)

    INPUT
    P1 = initial pressure (Pa)
    T1 = initial temperature (K)
    q = string of reactant species mole fractions
    mech = cti file containing mechanism data (i.e. 'gri30.cti')

    OUTPUT
    cj_speed = CJ detonation speed (m/s)
    R2 = R-squared value of LSQ curve fit

    """
    # DECLARATIONS
    numsteps = 20
    maxv = 2.0
    minv = 1.5

    w1 = np.zeros(numsteps+1, float)
    rr = np.zeros(numsteps+1, float)

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

    # INITIALIZE ERROR VALUES & CHANGE VALUES
    ERRFT = 1e-4
    ERRFV = 1e-4

    counter = 1
    r_squared = 0.0
    cj_speed = 0.0
    dnew = 0.0

    def curve_fit_function(x, a, b, c):
        return a * x**2 + b * x + c

    while (counter <= 4) or (r_squared < 0.99999):
        step = (maxv-minv)/float(numsteps)
        i = 0
        x = minv
        while x <= maxv:
            working_gas.TPX = [
                initial_temperature,
                initial_pressure,
                species_mole_fractions
                ]
            [working_gas,
             temp] = state_and_speed(
                 working_gas, initial_state_gas,
                 ERRFT,
                 ERRFV,
                 x
                 )
            w1[i] = temp
            rr[i] = working_gas.density / initial_state_gas.density
            i += 1
            x += step
        # [a, b, c, R2, SSE, SST] = LSQ_CJspeed(rr, w1)
        # Get curve fit and R^2
        [popt, _] = curve_fit(curve_fit_function, rr, w1)
        residuals = w1 - curve_fit_function(rr, *popt)
        ss_residual = np.sum(residuals**2)
        ss_total = np.sum((w1 - np.mean(w1))**2)
        r_squared = 1 - (ss_residual / ss_total)

        dnew = -popt[1] / (2. * popt[0])
        minv = dnew - dnew*0.001
        maxv = dnew + dnew*0.001
        counter += 1

    cj_speed = curve_fit_function(dnew, *popt)
    return [cj_speed, r_squared]
