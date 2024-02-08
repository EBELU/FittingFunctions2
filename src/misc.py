import numpy as np

def is_iter(obj):
    """
    Returns true if an object is iterable, else false.
    """
    try:
        iter(obj)
        return True
    except(TypeError):
        return False

def slice_spect(spect, *args, low = None, high = None):
    spect = np.array(spect)
    if (low is not None) and (high is not None):
        region = (spect >= low) & (spect <= high)
    elif (low is not None):
        region = (spect >= low)
    elif (high is not None):
        region = (spect <= high)
    else:
        raise RuntimeError("No limits given")
        return None
    
    if not args:
        return spect[region]
    else:
        spectra = [spect[region]]
        for spectrum in args:
            spectrum = np.array(spectrum)
            spectra.append(spectrum[region])
        return spectra
