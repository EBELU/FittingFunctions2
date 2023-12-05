#!/usr/bin/env python3
import numpy as np                       ## for handling of data
import matplotlib.pyplot as plt
import warnings
from scipy.optimize import curve_fit, OptimizeWarning


def GaussFunc(x, A, mu, sigma):
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))

def LineFunc(x, k, m):
    return k*x+m

class Peak:
    """A class to hold information regarding a single peak fitted using a Gussian"""

    def __init__(self, x, y,  region_start, region_stop):

        # these define a region of interest in our data binning:
        region = (region_start < x) & (x < region_stop)

        # now limit our data by 'slicing' it, limiting it to the region of interest:
        self.X = x[region]
        self.Y = y[region]

        A_guess = self.Y.max()
        mu_guess = self.X[np.where(self.Y == A_guess)[0][0]]
        #If you do the real FWHM you will crash if your peaks are too narrow
        FWHM = self.X[self.Y > A_guess/4]
        sigma_guess = FWHM[-1] - FWHM[0]

        guess = [A_guess, mu_guess, sigma_guess] # our initial guess of parameters for a Gaussian fit 

        ## scypi gives a warning if the fit does not work; we want to know about those, so we set them up to be caught here:
        with warnings.catch_warnings():
            warnings.simplefilter("error", OptimizeWarning)
            ## perform the gaussian fit to the data:
            try:
                ## use the scipy curve_fit routine (uses non-linear least squares to perform the fit)
                ## see http://docs.scipy.org/doc/scipy-0.16.1/reference/generated/scipy.optimize.curve_fit.html
                (self.A, self.mu, self.sigma), self.covar_matrix = curve_fit(GaussFunc,
                                        self.X,
                                        self.Y,
                                        p0=guess)
            except (RuntimeError, OptimizeWarning, TypeError):
                print("Gaussian fit failed! Try specifying another region.")

    
    def plot(self):
        #Choose colour of the Gaussian depending on if fit was okay
        color = 'forestgreen' if (self.covar_matrix[2][2]<self.sigma) else 'r'

        plt.plot(self.X, GaussFunc(self.X, self.A, self.mu, self.sigma), color=color, label = 'Gaussian fit')
        plt.legend(loc='upper right', frameon=False)


    def print_full_info(self):
        print("Estimated parameters:\n A = {:.5f}, mu = {:.5f},  sigma = {:.5f} \n".format(self.A, self.mu, self.sigma))
        print("Uncertainties in the estimated parameters: \n \u03C3\u00b2(A) = {:.5f}, \u03C3\u00b2(mu) = {:.5f}, \u03C3\u00b2(sigma) = {:.5f} \n".format(self.covar_matrix[0][0], self.covar_matrix[1][1], self.covar_matrix[2][2]))
        print("Covariance matrix: \n {}".format(self.covar_matrix))
        print("\n ======================================================= \n")


def perform_Gaussian_fit(x, y, regions, plotting = True, printing = True):
    peaks = []
    for region in regions:
        peaks.append(Peak(x, y, region[0], region[1]))

    if plotting:
        plt.figure()
        plt.step(x, y, where='mid', color='cornflowerblue', label='data') #plotting data

        for peak in peaks:
            peak.plot()

        plt.show()

    if printing:
        for peak in peaks:
            peak.print_full_info()

