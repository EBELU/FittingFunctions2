import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.signal import peak_widths

def GaussFunc(x, A, mu, sigma):
    return A*np.exp(-(x-mu)**2/(2*sigma**2))

def DGaussFunc(x, A1, mu1, sigma1, A2, mu2, sigma2):
    return GaussFunc(x, A1, mu1, sigma1) + GaussFunc(x, A2, mu2, sigma2)

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
    
    # The point closest to the peak in the gradient passing the threshold is
    # returend as the width.
    left_point = x_region[left_region][(left_grad <= threshold)][-1]
    right_point = x_region[right_region][(right_grad >= -threshold)][0]
    return (left_point, right_point)
        
"""Class to hold data for the gaussian"""
class gaussian:
    def __init__(self, A, mu, sigma, covar_matrix, x, y, region_limits=[None, None],
                 corr_f=lambda x: x * 0, corr_points=[None, None], areas=[0, 0, 0]):
        self.A, self.mu, self.sigma = A, mu, sigma
        self.covar_matrix = covar_matrix
        self._x, self._y = x, y
        self._region_limits = region_limits
        self._corr_f, self._corr_points = corr_f, corr_points
        self.G, self.N, self.B = areas

    "Getter functions"
    def value(self, x):
        return GaussFunc(x, self.A, self.mu, self.sigma)

    def area(self): #Analytical area calculation
        return self.A * abs(self.sigma) * np.sqrt(2 * np.pi)

    def FWHM(self, uncertanty=False):
        if uncertanty:
            return 2.35 * np.array([self.sigma, np.sqrt(self.covar_matrix[2][2])])
        else:
            return 2.35 * self.sigma

    def max_height(self):
        return self.A + self._corr_f(self.mu)


    "Plotting"
    def plot(self, xlabel=None, ylabel=None):
        plt.figure("peak", figsize=[9, 6])
        plt.plot(self._x, self._y, label="Data", color="deepskyblue")
        plt.plot(self._x, self._corr_f(self._x) + self.value(self._x),
                 label="Fitted Gaussian", color="indigo")
        if self._corr_points != [None, None]:
            try:
                plt.plot(self._x, self._corr_f(self._x), ":",
                         label="Scatter correction", color="red")
                plt.vlines(self._corr_points[0], 0, (self._corr_f(
                    self.mu)+self.A) * 1, color="red")
                plt.vlines(self._corr_points[1], 0, (self._corr_f(
                    self.mu)+self.A) * 1, color="red")
            except(AttributeError):
                pass
        plt.xlabel(xlabel, fontsize="large")
        plt.ylabel(ylabel, fontsize="large")
        plt.legend(fontsize="large")
        plt.xticks(fontsize="large")
        plt.yticks(fontsize="large")
        plt.show()

    def __str__(self):
        line = "========================= Gaussian Peak =========================\n"
        est_par = "Estimated parameters: A = {}, mu = {}, sigma = {}".format(
            round(self.A, 4), round(self.mu, 4), round(self.sigma, 4))
        uncert = "Uncertainties: \u03C3(A) = {}, \u03C3(mu) = {}, \u03C3(sigma) = {}".format(
            round(np.sqrt(self.covar_matrix[0][0]), 4),
            round(np.sqrt(self.covar_matrix[1][1]), 4),
            round(np.sqrt(self.covar_matrix[2][2]), 4))
        add_info = "Additional info: Max height = {}, G = {}, N = {}, B = {}".format(
            round(self.A + self._corr_f(self.mu), 4), round(self.G),
            round(self.N), round(self.B))
        covarmat = "Covariance matrix: \n {}".format(self.covar_matrix)
        return line + est_par + "\n" + uncert + "\n" + add_info + "\n\n" + covarmat + "\n\n"


def fit_gaussian(X, Y, region_start, region_stop,  # Setting region with one peak of interest
                 corr_left=None, corr_right=None,  # Scatter correction region
                 mu_guess=None, A_guess=None, sigma_guess=None,  # Manual set guesses
                 scatter_corr="auto", scatter_corr_points=3, corr_thresh=0.05):  # Scatter correction setting
    if corr_left != None and corr_right != None:
        # Warning for incorrectly placed correction limits
        if corr_left < region_start or corr_right > region_stop:
            print('\033[0;1;31m' + "Warning!" + '\033[0m')
            print("Scatter point outside region!\n")

    # Slice up region of interest
    X, Y = np.array(X), np.array(Y)
    region = (region_start < X) & (X < region_stop)
    x_region = X[region]
    y_region = Y[region]
    const_y_region = y_region.copy()

    # Make  correction for background scatter
    if scatter_corr:
        if scatter_corr == "auto":
            try:
                # Try to get boundries for correction, the funciton can be finiky
                corr_bounds = width_custom(
                    y_region, x_region, corr_thresh)
                if not corr_left:
                    corr_left = corr_bounds[0]
                if not corr_right:
                    corr_right = corr_bounds[1]
            except:
                corr_left, corr_right = (region_start, region_stop)
                print('\033[0;1;31m' + "Warning!" + '\033[0m')
                print(
                    f"Automatic scatter correction failed on the region [{region_start}, {region_stop}]")
                print("Try selecting another region or changing the threshold\n")
        else:
            if not corr_left:
                corr_left = region_start
            if not corr_right:
                corr_right = region_stop
                
        # The scatter_corr_points nr of data points inside the boundries for
        # scatter correction and average them in both x- and y-axis.
        left_xselection = sum(
            X[X >= corr_left][:scatter_corr_points])/scatter_corr_points
        right_xselection = sum(
            X[X <= corr_right][-scatter_corr_points:])/scatter_corr_points
        left_yselection = sum(
            Y[X >= corr_left][:scatter_corr_points])/scatter_corr_points
        right_yselection = sum(
            Y[X <= corr_right][-scatter_corr_points:])/scatter_corr_points

        # Make a linear fit to
        k, m = np.polyfit([left_xselection, right_xselection], [
                          left_yselection, right_yselection], 1)

        # Make a linear function from the parameters
        def corr_f(x): return k * x + m
        # Subtract the background region from the y-region
        y_region = const_y_region.copy() - corr_f(x_region)

        lin_region = ((x_region >= left_xselection) &
                      (x_region <= right_xselection))
        # Save correction data for plotting purposes
        Areas = [sum(const_y_region[lin_region]), sum(
            y_region[lin_region]), sum(corr_f(x_region[lin_region]))]
        corr_points = [corr_right, corr_left]

    if scatter_corr == False:
        def corr_f(x): return x * 0
        corr_points = [None, None]
        Areas = [sum(y_region), sum(y_region), 0.]

    # Assign guesses
    if not A_guess:
        # The guess for A is automatically set to the largest y-value in the region.
        A_guess = max(y_region)

    if not mu_guess:
        # The guess for mu is automatically set to the value on the x-axis corresponsing
        # to the largest y-value in the region.
        mu_guess = x_region[y_region.argmax()]

    if not sigma_guess:
        # The guess for sigma is set automatically to 1/2 of the width of the peak
        # at half max.
        sigma_guess = peak_widths(y_region, [y_region.argmax()], 0.5)[0][0]

    guesses = [A_guess, mu_guess, sigma_guess]

    # Handling the warning from SciPy as errors
    import warnings
    from scipy.optimize import OptimizeWarning
    with warnings.catch_warnings():
        warnings.simplefilter("error", OptimizeWarning)
        # perform the gaussian fit to the data:
        try:
            estimates, covar_matrix = curve_fit(GaussFunc,
                                                x_region,
                                                y_region,
                                                guesses)

            # create a Gauss object with the fitted coefficients for better code readability
            return gaussian(estimates[0], estimates[1], abs(estimates[2]), covar_matrix, x_region, const_y_region,
                            [region_start, region_stop], corr_f, corr_points, Areas)
        except (RuntimeError, OptimizeWarning, TypeError):
            raise RuntimeError(
                f"Gaussian fit failed at [{region_start}, {region_stop}]! Try specifying another region.")


def calibrate(Y, peak_regions, energies, plot=False, gauss=True):
    if len(peak_regions) != len(energies):
        raise IndexError("Nr of peaks must match number of energies provided!")
    X = np.arange(0, len(Y))
    peaks = []
    gauss_peaks = []
    for i, p_region in enumerate(peak_regions):
        if type(p_region).__name__ == "gaussian":
            gauss_peaks.append(p_region)
            peaks.append(p_region.mu)
            peak_regions[i] = p_region._region_limits
        elif gauss:
            gauss_peak = fit_gaussian(
                X, Y, region_start=p_region[0], region_stop=p_region[1], scatter_corr=True)
            gauss_peaks.append(gauss_peak)
            peaks.append(gauss_peak.mu)
        else:
            region = (p_region[0] < X) & (X < p_region[1])
            x_region = X[region]
            y_region = Y[region]
            peaks.append(x_region[y_region.argmax()])

    k, m = np.polyfit(peaks, energies, 1)
    def lin_f(x): return k * x + m
    calib_x = lin_f(X)
    peaks = np.array(peaks)

    if plot:
        plt.figure(figsize=[9, 6])
        plt.plot(calib_x, Y, color="lightskyblue", label="Data")
        for i, (start, stop) in enumerate(peak_regions):
            y_val = Y[round(peaks[i])]
            plt.vlines(calib_x[start], 0, y_val,
                       color="red", label="Region Limits")
            plt.vlines(calib_x[stop], 0, y_val, color="red")
            if gauss:
                region = (start < X) & (X < stop)
                gauss_x = calib_x[region]
                gauss_y = gauss_peaks[i].value(
                    X[region]) + gauss_peaks[i]._corr_f(X[region])
                plt.plot(gauss_x, gauss_peaks[i]._corr_f(
                    X[region]), color="red", label="Scatter correction", linestyle=":")
                plt.plot(gauss_x, gauss_y, color="indigo",
                         label="Fitted gaussian")
                plt.scatter(lin_f(gauss_peaks[i].mu), gauss_peaks[i].max_height(), color="teal", label="Calibration peak")
            else:
                plt.scatter(lin_f(peaks[i]), y_val,
                            color="teal", label="Calibration peak")
        plt.legend(fontsize="large")
        plt.xticks(fontsize="large")
        plt.yticks(fontsize="large")
        handles, labels = plt.gca().get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        plt.legend(by_label.values(), by_label.keys(), fontsize="large")
        plt.show()
    return (calib_x, (k, m))

def fit_double_gaussian(X, Y, region_start, region_stop, split_point, # Setting region with one peak of interest
                  corr_left=None, corr_right=None,  # Scatter correction region
                  plotting = False,
                  mu1_guess=None, mu2_guess=None, # Manual set guesses
                  A1_guess=None, A2_guess=None,
                  sigma1_guess=None, sigma2_guess=None,  
                  scatter_corr=True, scatter_corr_points=3):  # Scatter correction setting
    if corr_left != None and corr_right != None:
        # Warning for incorrectly placed correction limits
        if corr_left < region_start or corr_right > region_stop:
            print('\033[0;1;31m' + "Warning!" + '\033[0m')
            print("Scatter point outside region!\n")

    # Slice up region of interest
    region = (region_start < X) & (X < region_stop)
    x_region = X[region]
    y_region = Y[region]
    const_y_region = y_region.copy()

    # Make  correction for background scatter
    if scatter_corr:
        if not corr_left:
            corr_left = region_start
        if not corr_right:
            corr_right = region_stop
        # The scatter_corr_points nr of data points inside the boundries for
        # scatter correction and average them in both x- and y-axis.
        left_xselection = sum(
            X[X >= corr_left][:scatter_corr_points])/scatter_corr_points
        right_xselection = sum(
            X[X <= corr_right][-scatter_corr_points:])/scatter_corr_points
        left_yselection = sum(
            Y[X >= corr_left][:scatter_corr_points])/scatter_corr_points
        right_yselection = sum(
            Y[X <= corr_right][-scatter_corr_points:])/scatter_corr_points

        # Make a linear fit to
        k, m = np.polyfit([left_xselection, right_xselection], [
                          left_yselection, right_yselection], 1)

        # Make a linear function from the parameters
        def corr_f(x): return k * x + m
        # Subtract the background region from the y-region
        y_region = const_y_region.copy() - corr_f(x_region)

    else:
        def corr_f(x): return x * 0

    # Assign guesses
    if not A1_guess:
        # The guess for A is automatically set to the largest y-value in the region.
        A1_guess = max(y_region[x_region <= split_point])
    if not A2_guess:
        # The guess for A is automatically set to the largest y-value in the region.
        A2_guess = max(y_region[x_region > split_point])
    if not mu1_guess:
        # The guess for mu is automatically set to the value on the x-axis corresponsing
        # to the largest y-value in the region.
        mu1_guess = x_region[y_region[x_region <= split_point].argmax()]
    if not mu2_guess:
        # The guess for mu is automatically set to the value on the x-axis corresponsing
        # to the largest y-value in the region.
        left = x_region[x_region > split_point]
        mu2_guess = left[y_region[x_region > split_point].argmax()]

    if not sigma1_guess:
        # The guess for sigma is set automatically to 1/2 of the width of the peak
        # at half max.
        sigma1_guess = peak_widths(y_region[x_region <= split_point], [
                                   y_region[x_region <= split_point].argmax()], 0.5)[0][0]
    if not sigma2_guess:
        # The guess for sigma is set automatically to 1/2 of the width of the peak
        # at half max.
        sigma2_guess = peak_widths(y_region[x_region > split_point], [
                                   y_region[x_region > split_point].argmax()], 0.5)[0][0]
    guesses = [A1_guess, mu1_guess, sigma1_guess,
               A2_guess, mu2_guess, sigma2_guess]
    
    # Handling the warning from SciPy as errors
    import warnings
    from scipy.optimize import OptimizeWarning
    with warnings.catch_warnings():
        warnings.simplefilter("error", OptimizeWarning)
        # perform the gaussian fit to the data:
        try:
            estimates, covar_matrix = curve_fit(DGaussFunc,
                                                x_region,
                                                y_region,
                                                guesses)

            # create a Gauss object with the fitted coefficients for better code readability
            fitted_peaks = (gaussian(estimates[0], estimates[1], estimates[2], covar_matrix[:3, :3], 
                             x_region, const_y_region, [region_start, region_stop], corr_f=corr_f),
                    gaussian(estimates[3], estimates[4], estimates[5], covar_matrix[-3:, -3:], 
                             x_region, const_y_region, [region_start, region_stop], corr_f=corr_f))
            for P in fitted_peaks:
                N = P.value(x_region)
                P.B = sum(P._corr_f(x_region)[N > 1])
                P.N = sum(N[N > 1])
                P.G = P.N + P.B
        except (RuntimeError, OptimizeWarning, TypeError):
            print(f"Gaussian fit failed at [{region_start}, {region_stop}]! Try specifying another region.")
            return 0
    if plotting:
        plt.figure(figsize=[9,6])
        plt.plot(x_region, const_y_region, color = "lightskyblue", label = "Data")
        plt.vlines([region_start, region_stop], 0, fitted_peaks[0].max_height(), color  = "red", label = "Region limits")
        plt.plot(x_region, corr_f(x_region) , linestyle = ":", color = "red", label = "Scatter correction")
        plt.plot(x_region, fitted_peaks[0].value(x_region) + fitted_peaks[0]._corr_f(x_region), color = "indigo", label = "Fitted Gaussian 1")
        plt.plot(x_region, fitted_peaks[1].value(x_region) + fitted_peaks[1]._corr_f(x_region), color = "forestgreen", label = "Fitted Gaussian 2")
        plt.plot(x_region, fitted_peaks[0].value(x_region) + fitted_peaks[1].value(x_region) + fitted_peaks[1]._corr_f(x_region), linestyle = "--", color = "orchid", label = "Peak sum")
        plt.legend(fontsize="large")
        plt.xticks(fontsize="large")
        plt.yticks(fontsize="large")
        handles, labels = plt.gca().get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        plt.legend(by_label.values(), by_label.keys(), fontsize="large")
        plt.show()
    return fitted_peaks

