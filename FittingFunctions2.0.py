#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec  3 03:11:19 2023

@author: Erik Ewald
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.signal import peak_widths

def GaussFunc(x, A, mu, sigma):
    return A*np.exp(-(x-mu)**2/(2*sigma**2))
    
class gaussian:
    """Class to hold data for the gaussian"""
    def __init__(self, A, mu, sigma, covar_matrix):
        self.A = A
        self.mu = mu
        self.sigma = sigma
        self.covar_matrix = covar_matrix
    
    
    
class spectrum:
    """Class to hold spectrum data and fitted functions"""
    def __init__(self, yData, xData = None):
        self.ydata = yData
        self.xdata = xData
        
    def fit_gaussian(self, region_start, region_stop, 
                     mu_guess = None, A_guess = None, sigma_guess = None):
        #Slice up region of interest
        region = (region_start < self.x) & (self.x < region_stop)
        x_region = self.x[region]
        y_region = self.y[region]
        
        if A_guess:
            A_g = A_guess
        else:
            A_guess = max(y_region)
            
        if mu_guess:
            mu_g = mu_guess
        else:
            mu_g = x_region[y_region == max(y_region)]
            
        if sigma_guess:
            sigma_g = sigma_guess
        else:
            sigma_g = peak_widths(y_region, mu_guess, 0.3)[0]/2
            
        guesses = [A_g, mu_g, sigma_g]
            
        import warnings
        from scipy.optimize import OptimizeWarning
        with warnings.catch_warnings():
            warnings.simplefilter("error", OptimizeWarning)
            ## perform the gaussian fit to the data:
            try:
                ## use the scipy curve_fit routine (uses non-linear least squares to perform the fit)
                ## see http://docs.scipy.org/doc/scipy-0.16.1/reference/generated/scipy.optimize.curve_fit.html
                estimates, covar_matrix = curve_fit(GaussFunc,
                                        x_region,
                                        y_region,
                                        p0=guesses)
                ## create a Gauss object with the fitted coefficients for better code readability
                fitted_gaussian = gaussian(estimates[0], estimates[1], abs(estimates[2]), covar_matrix)
                return fitted_gaussian
            except (RuntimeError, OptimizeWarning, TypeError):
                print("Gaussian fit failed! Try specifying another mu_guess.")
                return 0
            
        
        
    