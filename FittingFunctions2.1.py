#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec  3 12:11:19 2023

@author: eewa
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.signal import peak_widths

def GaussFunc(x, A, mu, sigma):
    return A*np.exp(-(x-mu)**2/(2*sigma**2))

def width_custom(y_region, x_region, threshold = 200):
        # Find the centre point value and index
        mu = x_region[y_region == max(y_region)][0]
        i_peak = y_region.argmax()
        
        # Scipy makes an approximation for the width
        width = peak_widths(y_region, [i_peak], rel_height = 0.5)[0][0]
        #Using the approximation two ragions are cut ot on each side of the peak
        right_region = x_region > mu + width
        left_region = x_region < mu - width
        
        # The gradient is calcualted
        right_grad = np.gradient(y_region[right_region], 2)
        left_grad = np.gradient(y_region[left_region], 2)
        
        # plt.plot(x_region[right_region], right_grad)
        # plt.plot(x_region[left_region], left_grad)
        
        # The point closest to the peak in the gradient passing the threshold is
        # returend as the width.
        left_point = x_region[left_region][(left_grad <= threshold)][-1]
        right_point = x_region[right_region][(right_grad >= -threshold)][0]
        return(left_point, right_point)
        
"""Class to hold data for the gaussian"""
class gaussian:
    """Constructor performs fit"""
    def __init__(self, X, Y, region_start, region_stop, #Setting region with one peak of interest
                     corr_left = None, corr_right = None, corr_thresh = 0.05, #Scatter correction region
                     mu_guess = None, A_guess = None, sigma_guess = None, #Manual set guesses
                     scatter_corr = "auto", scatter_corr_points = 3):# Scatter correction setting
    
        if corr_left < region_start or corr_right > region_stop: # Warning for incorrectly placed correction limits
            print('\033[0;1;31m' + "Warning!" +'\033[0m')
            print("Scatter point outside region!")
            
        #Slice up region of interest
        region = (region_start < X) & (X < region_stop)
        x_region = X[region]
        y_region = Y[region]
        const_y_region = y_region
        
        # Make  correction for background scatter
        if scatter_corr:
            if scatter_corr == "auto" and not corr_left and not corr_right:
                try:
                    #Try to get boundries for correction, the funciton can be finiky
                    corr_left, corr_right = width_custom(y_region, x_region, corr_thresh)
                except:
                    corr_left, corr_right = (0,0)
                    print('\033[0;1;31m' + "Warning!" +'\033[0m')
                    print(f"Automatic scatter correction failed on the region [{region_start}, {region_stop}]")
                    print("Try selecting another region or changing the threshold")
            else:
                if not corr_left:
                        corr_left = region_start
                if not corr_right:
                        corr_right = region_stop
                    
            # The scatter_corr_points nr of data points outside the boundries for
            # scatter correction and average them in both x- and y-axis.
            left_xselection = sum(X[X <= corr_left][-scatter_corr_points:])/scatter_corr_points
            right_xselection = sum(X[X >= corr_right][:scatter_corr_points])/scatter_corr_points
            left_yselection = sum(Y[X <= corr_left][-scatter_corr_points:])/scatter_corr_points
            right_yselection = sum(Y[X >= corr_right][:scatter_corr_points])/scatter_corr_points
            
            # Make a linear fit to
            k, m = np.polyfit([left_xselection, right_xselection], [left_yselection, right_yselection], 1)
            
            # print([left_xselection, right_xselection])
            # print([left_yselection, right_yselection])
            # print(k, m)
            
            # Make a linear function from the parameters
            corr_f = lambda x : k * x + m
            # Subtract the background region from the y-region
            y_region = const_y_region.copy() - corr_f(x_region)
            
            lin_region = ((x_region >= left_xselection) & (x_region <= right_xselection))
            # Save correction data for plotting purposes
            self._corr_f = corr_f
            self._corr_points = [corr_left, corr_right]
            self._areas = [sum(const_y_region[lin_region]),sum(y_region[lin_region]), sum(self._corr_f(x_region[lin_region]))]
        if scatter_corr == False:
            self._corr_f = lambda x : x * 0
            self._corr_points = [None, None]
            self._areas = [sum(y_region),sum(y_region), sum(self._corr_f(x_region))]
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

        # plt.figure()
        # width = peak_widths(y_region, [y_region.argmax()], rel_height = 1)
        # print(width[1])
        # width = width[0][0]*2
        # test_1 = x_region > mu_g + width
        # test_2 = x_region < mu_g - width
        # plt.plot(x_region, ((y_region)))
        # plt.plot(x_region, (np.gradient(y_region, 2)))
        # plt.hlines(0, 1050, 1250)
        # plt.plot(x_region[y_region.argmax()], y_region[y_region.argmax()], ".")
        # plt.plot(x_region[test_1], y_region[test_1])
        # plt.plot(x_region[test_2], y_region[test_2])
        

        # plt.plot(x_region, y_region)
        
        # print(y_region.argmax())
        # print(sigma_g)
        # print(y_region[8])
        # print(guesses)
        
        # Handling the warning from SciPy as errors
        import warnings
        from scipy.optimize import OptimizeWarning
        with warnings.catch_warnings():
            warnings.simplefilter("error", OptimizeWarning)
            ## perform the gaussian fit to the data:
            try:
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
                self._region = [x_region, const_y_region]
                self._region_limits = [region_start, region_stop]
            except (RuntimeError, OptimizeWarning, TypeError):
                raise RuntimeError(f"Gaussian fit failed at [{region_start}, {region_stop}]! Try specifying another region.")

                
    "Getter functions"
    def value(self, x):
        return GaussFunc(x, self.A, self.mu, self.sigma)
    def area(self):
        return self.A*np.abs(self.sigma)*np.sqrt(2*np.pi)
    def FWHM(self):
        return 2.35 * np.array([self.sigma, np.sqrt(self.covar_matrix[2][2])])
    def addinfo(self):
        return [self.A + self._corr_f(self.mu), self._areas[0], self._areas[1],self._areas[2]]
    
    "Plotting"
    def plot(self, xlabel = None, ylabel = None):
        plt.figure(figsize = [9,6])
        plt.plot(self._region[0], self._region[1], label = "Data", color = "deepskyblue")
        plt.plot(self._region[0], self._corr_f(self._region[0]) + self.value(self._region[0]), label = "Fitted Gaussian", color = "forestgreen")
        if self._corr_points != [None, None]:
            try:
                plt.plot(self._region[0], self._corr_f(self._region[0]), label = "Scatter correction", color = "crimson")
                plt.vlines(self._corr_points[0], 0, (self._corr_f(self.mu)+self.A) * 1, color = "crimson")
                plt.vlines(self._corr_points[1], 0, (self._corr_f(self.mu)+self.A) * 1, color = "crimson")
            except(AttributeError):
                pass
        plt.xlabel(xlabel,fontsize="large")
        plt.ylabel(ylabel,fontsize="large")
        plt.legend(fontsize="large")
        plt.xticks(fontsize="large")
        plt.yticks(fontsize="large")

    def __str__(self):
        est_par = "Estimated paramters: A = {}, mu = {}, sigma = {}".format(
             round(self.A, 4), round(self.mu, 4), round(self.sigma, 4))
        uncert = "Uncertanties: \u03C3(A) = {}, \u03C3(mu) = {}, \u03C3(sigma) = {}".format(
            round(np.sqrt(self.covar_matrix[0][0]), 4), 
            round(np.sqrt(self.covar_matrix[1][1]), 4), 
            round(np.sqrt(self.covar_matrix[2][2]), 4))
        add_info = "Addtional info: Max height = {}, G = {}, N = {}, B = {}".format(
            round(self.A + self._corr_f(self.mu), 4), round(self._areas[0], 4),
            round(self._areas[1], 4), round(self._areas[2], 4))
        covarmat = "Covariance matrix: \n {}".format(self.covar_matrix)
        return est_par + "\n" + uncert + "\n" + add_info + "\n\n" + covarmat + "\n\n"






