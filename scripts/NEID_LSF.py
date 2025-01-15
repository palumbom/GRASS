from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from scipy.special import erf

import sys
import os

# Add the directory containing the script to sys.path
sys.path.append(os.path.abspath(os.path.dirname(__file__)))

def tophat_gauss_kernel( x, mean, fwhm, box_width ):
    """A function to return a top hat convolved Gaussian kernel provided the two width parameters.
    Note that the kernel values are not normalized (needed for convolution).
    
    Parameters
    ----------
    x : array
        An array of pixel locations to compute the kernel at.
    
    mean : float
        The mean of the kernel at the point in the observed spectrum to convolve at.
        
    fwhm : float
        The Gaussian FWHM of the kernel at the point in the observed spectrum to convolve at.
        
    box_width : float
        The top hat box width of the kernel at the point in the observed spectrum to convolve at.
    
    Returns
    -------
    kernel_values : array
        The values of the LSF convolution kernel.
    """
    
    # First convert the FWHM to a sigma
    sigma = fwhm / ( 2 * np.sqrt( 2 * np.log( 2 ) ) )
    
    # Generate the arguments for the two erf components.
    erf_arg_1 = ( ( 2 * mean + box_width ) - 2 * x ) / ( 2 * np.sqrt( 2 ) * sigma )
    erf_arg_2 = ( (-2 * mean + box_width ) + 2 * x ) / ( 2 * np.sqrt( 2 ) * sigma )
    
    # Generate the kernel values
    kernel_values = erf( erf_arg_1 ) + erf( erf_arg_2 )
    
    return kernel_values

def read_lsf_masterfile(file_name):
    """Reads in a text file (in form of IRAF beam trace file) as a dictionary.
    
    Taken from the PRVextract.libs.beamtrace library file in the NEID DRP.

    Parameters
    ----------
    file_name : str
        The path to the text file to read in.

    Returns
    -------
    file_dict : dict of dicts
        The dictionary corresponding to the input file.
    """
    
    # Initialize output dictionary
    file_dict = {}

    # Go through each of the lines in the file to append to output dictionary
    current_aperture = None
    for line in open(file_name):
        line = line.split()

        # Skip empty lines
        if not line:
            continue

        # Read the line
        key = line[0]
        value = [parseValue(value) for value in line[1:]]

        # If only one value for line, turn from list to single value
        if len(value) == 1:
            value = value[0]

        # Insert the information into the dictionary
        if key == 'APERTURE':
            if current_aperture is not None:
                file_dict[current_aperture['APERTURE']] = current_aperture
            current_aperture = dict({'APERTURE': value})
        else:
            current_aperture[key] = value

    if current_aperture is not None:
        file_dict[current_aperture['APERTURE']] = current_aperture

    return file_dict

def parseValue(value):
    """Utility function for parsing a value when reading in beam trace file.
    """
    
    for parse_function in [int, float, str]:
        try:
            return parse_function(value)
        except:
            continue
            
    return value

def get_variable_lsf_kernel_values(obs_pixel, kernel_pixel_arr, wavelength_spacing, order):
    """Function to generate a 1D LSF kernel array given parameterization defined in the input masterfile.
    The LSF profile has parameters that are fit with polynomials in pixel space across each order.
    
    Parameters
    ----------
    obs_pixel : float
        The pixel location within spectral order to generate LSF at.
    
    kernel_pixel_arr : array.
        The pixel locations at which the kernel will be evaluated.
    
    kernel_profile_information : dict
        Dictionary with GLOBAL information about the LSF profile parameterization.
    
    kernel_parameters : dict
        The order-specific dictionary from the input LSF masterfile containing the kernel parameter fit information.
    
    wavelength_spacing : float
        The equal spacing of the wavelength grid being convolved, for normalizing the kernel.
    
    Returns
    -------
    lsf_kernel_values : array
        The values of the LSF convolution kernel for the input kernel_pixel_arr.
    """
    lsf_info = read_lsf_masterfile("data/neidMaster_LSFParameterCoefficients20230824HRSci_v1.txt")
    kernel_profile_information = lsf_info['GLOBAL']
    kernel_parameters = lsf_info[order]
    # Scale the input observation pixel value to be in the domain over which the LSF parameters are fit (turning [0,9215] to [-1,1])
    obs_pixel_m1top1 = scale_interval_m1top1(obs_pixel, *kernel_profile_information['DOMAIN'])
    
    ### Pick out the right line profile parameterization from the input master file. Currently just top hat gaussian
    if kernel_profile_information['PROFILE'] == 'TOP_HAT_GAUSSIAN':
        
        # Take the polynomial coefficients describing the kernel parameters vs. pixel for this order
        kernel_fwhm_poly = kernel_parameters['COEFF_1'] # Coefficient 1 is the Gaussian FWHM
        kernel_boxwidth_poly = kernel_parameters['COEFF_2'] # Coefficient 2 is the top hat box width
        
        # Now apply the fits to the input observed pixel value
        kernel_fwhm = np.polyval( kernel_fwhm_poly, obs_pixel_m1top1 )
        kernel_boxwidth = np.polyval( kernel_boxwidth_poly, obs_pixel_m1top1 )
        
        # Generate the actual kernel values
        lsf_kernel_values = tophat_gauss_kernel(kernel_pixel_arr, obs_pixel, kernel_fwhm, kernel_boxwidth)
        
        # Normalize the kernel values so that their integral is 1 (for flux conservation)
        lsf_kernel_values /= np.trapz(lsf_kernel_values, dx = wavelength_spacing)
        
    return lsf_kernel_values

def scale_interval_m1top1(x,a,b):
    """Scales input x values over interval a to b onto the range of -1 to 1.
    """
    
    return (2.0 * x - (b + a)) / (b - a)
