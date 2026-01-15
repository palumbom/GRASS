import numpy  as np
import numpy as np
import matplotlib.pyplot as plt
import pickle

class NeidLSFModel:
    """
    This class models the NIED Instruemt Line Spread
    Function Profile.

    Only files needed are this file and "NEID_LSFMODEL_##_###".
    
    Parameters:
    -----------
    lsf_params_dict : Content of the pickle file "lsf_params_pickle".

    Notes:
    ------
        -Run the function "lsf_model(echelle_order, pixel_idx)" to get LSF.
        -"echelle_order" : Always refers to python index of a row in the Neid Spectrum.
        -This implementation spans the "echelle_orders" from order_start to order_end.
        -This implementation is expected to behave best within the Free Spectral Range.
    """

    def __init__(self, lsf_params_dict):
        self.lsf_params_dict = lsf_params_dict
        self.sciwave = self.lsf_params_dict["sciwave"]
        self.order_start = self.lsf_params_dict["order_start"]
        self.order_end = self.lsf_params_dict["order_end"]
        self.lam_ref = self.lsf_params_dict["lam_ref"]
        self.width_ref = self.lsf_params_dict["width_ref"]
        self.echelle_order_ref = self.lsf_params_dict["echelle_order_ref"]
        self.pixel_idx_ref = self.lsf_params_dict["pixel_idx_ref"]
        self.sigma_gaussian = self.lsf_params_dict["sigma_gaussian"]
        self.ortho_basis = self.lsf_params_dict["ortho_basis"]
        self.legendre_basis = self.lsf_params_dict["legendre_basis"]
        self.legendre_coeffs = self.lsf_params_dict["legendre_coeffs"]
        self.grad_sciwave = np.gradient(self.sciwave, axis=1)
  
            
    def convolve(self, x1, y1, x2, y2):
        """
        This function convolves two kernels in "full" mode.
        
        Parameters:
        -----------
        x1 : array
            x points for signal.
        y1 : array
            y points for signal.
        x2 : array
            x points for kernel.
        y2 : array
            y points for kernel.

        Returns:
        --------
        y_conv : array
            Convolved signal sampled over x1.
        """

        dx1 = x1[1]-x1[0]
        dx2 = x2[1]-x2[0]
        
        tol=1e-6

        if np.abs(dx1-dx2)>tol:
            return "Sampling rate is different in both kernels."
        
        y_conv = np.convolve(y1, y2)
        x_conv = np.linspace(x1[0]-np.abs(x2[0]), x1[-1]+np.abs(x2[-1]), len(y_conv))

        y_conv = np.interp(x1, x_conv, y_conv)

        return y_conv
        
    
    def gaussian(self, pixel):
        """
        Gaussian spread function representing the "detector-LSF".
        Note that the kernel should be normalized before use.

        Parameters:
        -----------
        pixel : float
            Pixel position.

        Returns:
        --------
        float :
            A gaussian profile as a function of pixel position.
        """
        return np.exp(-(pixel**2)/(2 * self.sigma_gaussian**2))
    

    def aperture_width(self, echelle_order, pixel):
        """
        Returns the width of the image of the optical fibre at 
        a given order and pixel.

        Parameters:
        -----------
        echelle_order : int
            Echelle Order.
        pixel : int
            Pixel index.

        Returns:
        --------
        float : 
            Width of image of optical fibre.
        """
        order = echelle_order - self.order_start
        lam = self.sciwave[order, pixel]
        lam_ref = self.lam_ref
        width_ref = self.width_ref
        grad_sciwave = self.grad_sciwave

        order_ref = self.echelle_order_ref - self.order_start
        pixel_ref = self.pixel_idx_ref

        width = (lam / lam_ref) * width_ref * \
            grad_sciwave[order_ref, pixel_ref] / grad_sciwave[order, pixel]
        return width


    def aperture_image(self, pixel, width):
        """
        The image of the circular aperture on the detector. The 
        width can be tuned using echelle order and pixel number.
        Note that the kernel should be normalized before use.
        
        Parameters:
        -----------
        pixel : float
            Pixel positions.
        width : float
            Width of the image in pixel units.
            
        Returns:
        float :
            Image profile (2 * semicircle).
        """

        if (np.abs(pixel)>=width):
            return 0
        y = 2 * np.sqrt(width**2 - pixel**2)

        return y
    

    def pixel_response(self, pixel):
        """
        A flat pixel response function. Note that the kernel 
        should be normalized before use.

        Parameters:
        -----------
        pixel : float
            Pixel positions.
        
        Returns:
        --------
        float :
            A uniform function with width = 1 pixel.
        """
        if (-0.5<=pixel<=0.5):
            return 1
        return 0


    def get_weight_function(self, echelle_order, pixel_index):
        """
        Returns the weight_function / "physical LSF" for a given order and pixel.
        
        Parameters:
        -----------
        echelle_order : int
            Echelle order
        pixel_index : int
            Pixel index
            
        Returns:
        --------
        pixel_vals : 
            Highly sampled pixel points.
        weights_model
            Weights sampled over pixel_vals.
        """

        width_aperture_image = self.aperture_width(echelle_order, pixel_index)
        pixel_vals = np.linspace(-width_aperture_image-4, width_aperture_image+4, 1001)

        gaussian_profile = self.gaussian(pixel_vals)
        gaussian_profile = gaussian_profile/np.sum(gaussian_profile)

        aperture_image_profile = np.array([self.aperture_image(pixel, width_aperture_image)
                                   for pixel in pixel_vals])
        aperture_image_profile = aperture_image_profile/np.sum(aperture_image_profile)
    
        pixel_response_profile = np.array([self.pixel_response(pixel) for pixel in pixel_vals])
        pixel_response_profile = pixel_response_profile / np.sum(pixel_response_profile)


        weights_model = self.convolve(pixel_vals, gaussian_profile, pixel_vals, aperture_image_profile)
        weights_model = self.convolve(pixel_vals, weights_model, pixel_vals, pixel_response_profile)

        weights_model = weights_model/np.sum(weights_model)

        return pixel_vals, weights_model
    

    def get_polycoeffs(self, echelle_order, pixel_idx):
        legendre_basis = self.legendre_basis
        legendre_coeffs = self.legendre_coeffs[echelle_order]
        leg_poly = [sum(c * P for c, P in zip(row, legendre_basis)) for row in legendre_coeffs]
        coeff_vals = np.array([lp(pixel_idx) for lp in leg_poly])
        return coeff_vals
    

    def lsf_model(self, echelle_order, pixel_idx):
        coeff_vals = self.get_polycoeffs(echelle_order, pixel_idx)
        pixels, weights = self.get_weight_function(echelle_order, pixel_idx)

        ortho_pvals = []
        for p in self.ortho_basis:
            ortho_pvals.append(p(pixels))

        model_vals = weights * (1 + np.sum(coeff_vals[:, None] * ortho_pvals, axis=0))
        model_vals /= np.sum(model_vals)

        return pixels, model_vals

def lsf_model_export(echelle_order, pixel_idx, wavelength, fluxes, ref_wavelength):
    pkl_filename = "data/NEID_LSFMODEL_HR_SCI"
    with open(pkl_filename, "rb") as f:
        lsf_params_dict = pickle.load(f)
    LSFMODEL = NeidLSFModel(lsf_params_dict)

    echelle_order = echelle_order - 1
    pixel_vals, lsf_vals = LSFMODEL.lsf_model(echelle_order, pixel_idx)

    sciwave_ref = lsf_params_dict["sciwave"]
    grad_sciwave_ref = np.gradient(sciwave_ref, axis=1)
    dlam_dpix = grad_sciwave_ref[echelle_order, pixel_idx]

    dpix = pixel_vals[1] - pixel_vals[0]
    lams = wavelength - ref_wavelength
    pixels = lams / dlam_dpix

    pixels_new  = np.arange(pixels[0], pixels[-1]+dpix, dpix)
    fluxes_new = np.interp(pixels_new, pixels, fluxes)

    fluxes_conv = 1 - np.convolve(1 - fluxes_new, lsf_vals, mode="same") 
    return pixels_new, fluxes_conv, dlam_dpix