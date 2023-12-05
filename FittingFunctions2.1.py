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

x = np.linspace(0, 100, 500)
y = GaussFunc(x, 5, 20, 3)

class gaussian:
    """Class to hold data for the gaussian"""
    def __init__(self, x, y, region_start, region_stop, #Setting region with one peak of interest
                     mu_guess = None, A_guess = None, sigma_guess = None, #Manual set guesses
                     scatter_corr = True, scatter_corr_points = 3, # Scatter correction setting
                     scatter_corr_left = None, scatter_corr_right = None):
        #Slice up region of interest
        region = (region_start < x) & (x < region_stop)
        x_region = x[region]
        y_region = y[region]
        
        #Assign guesses
        if A_guess:
            A_g = A_guess
        else:
            #The guess for A is automatically set to the largest y-value in the region.
            A_g = max(y_region)
            
        if mu_guess:
            mu_g = mu_guess
        else:
            # The guess for mu is automatically set to the value on the x-axis corresponsing
            # to the largest y-value in the region.
            mu_g = x_region[y_region == max(y_region)][0]

        if sigma_guess:
            sigma_g = sigma_guess
        else:
            # The guess for sigma is set automatically to 1/2 of the width of the peak
            # at half max.
            sigma_g = peak_widths(y_region, [y_region.argmax()], 0.5)[0][0]
            
        guesses = [A_g, mu_g, sigma_g]
        
        """MÖG FRÅN DEBUGGING"""
        #plt.plot(x_region, y_region)
        
        # print(y_region.argmax())
        # print(sigma_g)
        # print(y_region[8])
        # print(guesses)
        
        #Fitting the parameters
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
                                        guesses)
                
               # print((estimates[0], estimates[1], abs(estimates[2])))
                ## create a Gauss object with the fitted coefficients for better code readability
                self.A = estimates[0]
                self.mu = estimates[1]
                self.sigma = abs(estimates[2])
                self.covar_matrix = covar_matrix
            except (RuntimeError, OptimizeWarning, TypeError):
                print("Gaussian fit failed! Try specifying another mu_guess.")

        
    def value(self, x):
        return GaussFunc(x, self.A, self.mu, self.sigma)
    def area(self):
        return self.A*np.abs(self.sigma)*np.sqrt(2*np.pi)

    def __str__(self):
        est_par = "Estimated paramters: A = {}, mu = {}, sigma = {}".format(
            round(self.A, 4), round(self.mu, 4), round(self.sigma, 4))
        uncert = "Uncertanties: \u03C3(A) = {}, \u03C3(mu) = {}, \u03C3(sigma) = {}".format(
            round(np.sqrt(self.covar_matrix[0][0]), 4), 
            round(np.sqrt(self.covar_matrix[1][1]), 4), 
            round(np.sqrt(self.covar_matrix[2][2]), 4))
        covarmat = "Covariance matrix: \n {}".format(self.covar_matrix)
        return est_par + "\n" + uncert +"\n\n" + covarmat + "\n\n"
