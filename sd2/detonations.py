# -*- coding: utf-8 -*-
"""
Functions for detonation calculations.

Original functions from Shock and Detonation Toolbox
http://www.galcit.caltech.edu/EDL/public/cantera/html/SD_Toolbox/
"""

import warnings
import numpy as np
import cantera as ct
import multiprocessing as mp
from scipy.optimize import curve_fit
from . import calculate_error, states


def calculate_cj_state(
        working_gas,
        initial_state_gas,
        error_tol_temperature,
        error_tol_velocity,
        density_ratio,
        max_iterations=500
):
    """
    This function calculates the Chapman-Jouguet state and wave speed using
    Reynolds' iterative method.

    Original function: CJ_calc in PostShock.py

    Parameters
    ----------
    working_gas : cantera.composite.Solution
        A cantera gas object used for calculations (???).
    initial_state_gas : cantera.composite.Solution
        A cantera gas object for the working gas mixture in its initial,
        undetonated state.
    error_tol_temperature : float
        Temperature error tolerance for iteration.
    error_tol_velocity : float
        Velocity error tolerance for iteration.
    density_ratio : float
        density ratio.
    max_iterations : int
        Maximum number of loop iterations used to calculate output. Default
        is 500.

    Returns
    -------
    working_gas : cantera.composite.Solution
        Gas object at equilibrium state.
    initial_velocity : float
        Initial velocity resulting in the input density ratio, in m/s.
    """
    # initial state
    initial_volume = 1 / initial_state_gas.density

    # set guess values
    guess_temperature = 2000
    guess_velocity = 2000
    guess_volume = initial_volume / density_ratio
    guess_density = 1 / guess_volume

    # set deltas
    delta_temperature = 1000
    delta_velocity = 1000

    # equilibrate
    states.get_equilibrium_properties(
        working_gas,
        guess_density,
        guess_temperature
    )

    loop_counter = 0
    while (
            abs(delta_temperature) > (error_tol_temperature * guess_temperature)
            or
            abs(delta_velocity) > (error_tol_velocity * guess_velocity)
    ):
        loop_counter += 1
        # check for non-convergence
        if loop_counter == max_iterations:
            warnings.warn(
                'No convergence within {0} iterations'.format(max_iterations)
            )
            return [working_gas, guess_velocity]

        # calculate unperturbed enthalpy and pressure error for current guess
        [error_enthalpy,
         error_pressure] = calculate_error.equilibrium(
            working_gas,
            initial_state_gas,
            guess_velocity
        )

        # perturb temperature
        delta_temperature = 0.02 * guess_temperature
        perturbed_temperature = guess_temperature + delta_temperature
        states.get_equilibrium_properties(
            working_gas,
            guess_density,
            perturbed_temperature
        )

        # calculate error rates for temperature perturbed state
        [error_perturbed_enthalpy,
         error_perturbed_pressure] = calculate_error.equilibrium(
            working_gas,
            initial_state_gas,
            guess_velocity
        )
        derivative_enthalpy_temperature = (error_perturbed_enthalpy -
                                           error_enthalpy) / delta_temperature
        derivative_pressure_temperature = (error_perturbed_pressure -
                                           error_pressure) / delta_temperature

        # perturb velocity
        delta_velocity = 0.02 * guess_velocity
        perturbed_velocity = guess_velocity + delta_velocity
        perturbed_temperature = guess_temperature
        states.get_equilibrium_properties(
            working_gas,
            guess_density,
            perturbed_temperature
        )

        # calculate error rates for velocity perturbed state
        [error_perturbed_enthalpy,
         error_perturbed_pressure] = calculate_error.equilibrium(
            working_gas,
            initial_state_gas,
            perturbed_velocity
        )
        derivative_enthalpy_velocity = (error_perturbed_enthalpy -
                                        error_enthalpy) / delta_velocity
        derivative_pressure_velocity = (error_perturbed_pressure -
                                        error_pressure) / delta_velocity

        # invert matrix
        jacobian = derivative_enthalpy_temperature *\
            derivative_pressure_velocity - \
            derivative_pressure_temperature * \
            derivative_enthalpy_velocity
        b = [derivative_pressure_velocity,
             -derivative_enthalpy_velocity,
             -derivative_pressure_temperature,
             derivative_enthalpy_temperature]
        a = [-error_enthalpy,
             -error_pressure]
        delta_temperature = (b[0] * a[0] + b[1] * a[1]) / jacobian
        delta_velocity = (b[2] * a[0] + b[3] * a[1]) / jacobian

        # limit temperature changes
        max_temperature_delta = 0.2 * guess_temperature
        if abs(delta_temperature) > max_temperature_delta:
            delta_temperature *= max_temperature_delta / abs(delta_temperature)

        # apply deltas and equilibrate
        guess_temperature += delta_temperature
        guess_velocity += delta_velocity
        states.get_equilibrium_properties(
            working_gas,
            guess_density,
            guess_temperature
        )

    return [working_gas, guess_velocity]


def _calculate_over_ratio_range(
        current_state_number,
        current_density_ratio,
        initial_temperature,
        initial_pressure,
        species_mole_fractions,
        mechanism,
        error_tol_temperature,
        error_tol_velocity
):
    initial_state_gas = ct.Solution(mechanism)
    initial_state_gas.TPX = [
        initial_temperature,
        initial_pressure,
        species_mole_fractions
    ]
    working_gas = ct.Solution(mechanism)
    working_gas.TPX = [
        initial_temperature,
        initial_pressure,
        species_mole_fractions
    ]
    [_,
     current_velocity] = calculate_cj_state(
        working_gas,
        initial_state_gas,
        error_tol_temperature,
        error_tol_velocity,
        current_density_ratio
    )

    return current_state_number, current_velocity


def calculate_cj_speed(
        initial_pressure,
        initial_temperature,
        species_mole_fractions,
        mechanism,
        use_multiprocessing=False,
        return_r_squared=False,
        return_state=False
):
    """
    This function calculates CJ detonation velocity

    Original function: CJspeed in PostShock.py

    Parameters
    ----------
    initial_pressure : float
        initial pressure (Pa)
    initial_temperature : float
        initial temperature (K)
    species_mole_fractions : str or dict
        string or dictionary of reactant species mole fractions
    mechanism : str
        cti file containing mechanism data (e.g. 'gri30.cti')
    use_multiprocessing : bool
        use multiprocessing to speed up CJ speed calculation
    return_r_squared : bool
        return the R^2 value of the CJ speed vs. density ratio fit
    return_state : bool
        return the CJ state corresponding to the calculated velocity

    Returns
    -------
    dict
    """
    # DECLARATIONS
    num_steps = 5
    max_density_ratio = 2.0
    min_density_ratio = 1.5

    if use_multiprocessing:
        pool = mp.Pool()

    # Set error tolerances for CJ state calculation
    error_tol_temperature = 1e-4
    error_tol_velocity = 1e-4

    counter = 1
    r_squared = 0.0
    delta_r_squared = 0.0
    adjusted_density_ratio = 0.0

    def curve_fit_function(x, a, b, c):
        """
        Quadratic function for least-squares curve fit of cj speed vs. density
        ratio
        """
        return a * x**2 + b * x + c

    while (counter <= 4) and (r_squared < 0.99999 or delta_r_squared < 1e-7):
        density_ratio_array = np.linspace(
            min_density_ratio,
            max_density_ratio,
            num_steps + 1
        )

        if use_multiprocessing:
            # parallel loop through density ratios
            stargs = [[number,
                       ratio,
                       initial_temperature,
                       initial_pressure,
                       species_mole_fractions,
                       mechanism,
                       error_tol_temperature,
                       error_tol_velocity
                       ]
                      for number, ratio in
                      zip(
                          range(len(density_ratio_array)),
                          density_ratio_array
                      )
                      ]
            result = pool.starmap(_calculate_over_ratio_range, stargs)

        else:
            # no multiprocessing, just use map
            result = list(map(
                _calculate_over_ratio_range,
                [item for item in range(len(density_ratio_array))],
                density_ratio_array,
                [initial_temperature for _ in density_ratio_array],
                [initial_pressure for _ in density_ratio_array],
                [species_mole_fractions for _ in density_ratio_array],
                [mechanism for _ in density_ratio_array],
                [error_tol_temperature for _ in density_ratio_array],
                [error_tol_velocity for _ in density_ratio_array]
            ))

        result.sort()
        cj_velocity_calculations = np.array(
            [item for (_, item) in result]
        )

        # Get curve fit
        [curve_fit_coefficients, _] = curve_fit(
            curve_fit_function,
            density_ratio_array,
            cj_velocity_calculations
            )

        # Calculate R^2 value
        residuals = cj_velocity_calculations - curve_fit_function(
            density_ratio_array,
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

    if return_state:
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
        calculate_cj_state(
            working_gas,
            initial_state_gas,
            error_tol_temperature,
            error_tol_velocity,
            adjusted_density_ratio
        )

    if return_r_squared and return_state:
        return {'cj speed': cj_speed,
                'R^2': r_squared,
                'cj state': working_gas
                }
    elif return_state:
        return {'cj speed': cj_speed,
                'cj state': working_gas
                }
    elif return_r_squared:
        return {'cj speed': cj_speed,
                'R^2': r_squared
                }
    else:
        return {'cj speed': cj_speed}
