# -*- coding: utf-8 -*-
"""
Created on Fri Apr 27 11:26:35 2018

@author: cartemic
"""


def perturb(property_str, property_dict, delta_dict):
    """
    This function perturps a property by some amount delta.

    Parameters
    ----------
    property_str : str
        A cantera gas object used for calculations (???).
    property_dict : dict
        A dictionary containing properties (as keys) and their numeric values
        (as values)
    delta_dict : dict
        A dictionary containing numeric property deltas with properties as keys

    Returns
    -------
    perturbed_dict : dict
        A dictionary containing properties as keys, with the appropriate value
        being perturbed by the corresponding delta
    """
    # Ensure correct input types

    # Ensure property_str is a key in both dictionaries

    perturbed_dict = property_dict.copy()
    perturbed_dict[property_str] += delta_dict[property_str]

    return perturbed_dict
