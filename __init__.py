"""
An improved version of fittingfunctions.

import FittingFunctions2 as ff

Contains functions for fitting Gaussian distributions in spectra.
Additional functions for calibration and slicing spectra are also included.

See accompanying README files for documentation on provided functions.
"""

from .src.gaussian_fitting import fit_gaussian, fit_double_gaussian, calibrate
from .src.misc import is_iter, slice_spect
from .src.UI import UI_wrapper as UI
