import numpy as np
from fit_cont._cont_src import fit_sdss_cont


## function to interpolate the output of fit_cont
def get_cont(C,wavearray):
    ressize = int(C[0])  ## extract params
    rr = np.arange(0,ressize,1)
    wavegood=C[rr+1]
    fluxgood=C[rr+1+ressize]
    y2=C[rr+1+2*ressize]
    continuum_array = np.zeros_like(wavearray)
    ib = (wavearray > wavegood[1]) & (wavearray < wavegood[-1])
    klo = np.zeros(wavearray.shape,dtype=int)
    khi = np.zeros(wavearray.shape,dtype=int)
    for j in np.arange(wavearray.shape[0])[ib]:  ## populate continuum array
        wavetest = wavearray[j]
        # cubic interpolation
        klo[j] = np.where(wavegood<wavetest)[0][-1]
        khi[j] = np.where(wavegood>wavetest)[0][0]
    h=wavegood[khi[ib]]-wavegood[klo[ib]];
    a=(wavegood[khi[ib]]-wavearray[ib])/h;
    b=(wavearray[ib]-wavegood[klo[ib]])/h;
    continuum_array[ib]=a*fluxgood[klo[ib]]+b*fluxgood[khi[ib]]+(((a*a*a)-a)*y2[klo[ib]]+((b*b*b)-b)*y2[khi[ib]])*(h*h)/6.
    return continuum_array


# fits a spline to spectrum; given wavelength array, flux array, flux errors, redshift
def spline_fit(wave,flux,sigma,z_in,dpx2):
    # Automated continuum fitter parameters. Built-in params for steepness and overall flux adjustment.
    # Avoid changing these.
    slopethresh = 0.033 
    fluxthresh = 0.99
    Fscale = 1.0

    # Normalisation wavelength
    wave_norm = 1290.0


    # Fit range
    wave_rest_array = np.arange(1000,3000,0.5)


    # Calculate rest-frame wavelength
    wave_rest = wave/(1+z_in)


    # Set up the arrays to hold the continuum
    cont =  np.zeros_like(flux)
    cont_norm = np.zeros_like(flux)
    flux_norm = np.zeros_like(flux)
    sigma_norm = np.zeros_like(flux)


    # First round of continuum fitting for a good first guess. Sufficient for lambda>1215 Ang.

    # Mask to avoid bad values
    ifit = wave > 0.0
    nfit = np.sum(ifit) # total number of fittable pixels MUST be provided to routine

    wave=wave[ifit]
    flux=flux[ifit]
    sigma=sigma[ifit]
    wave_rest=wave_rest[ifit]


    # Most important fitting parameters.
    # QSOs with narrow lines benefit from decreasing deltapix+minpix
    # But beware overfitting of broad absorption features!
    # deltapix1 is the spline window redward of Lya, while deltapix2 is blueward.
    deltapix1 = 5   ## OG parameters for first round
    deltapix2 = dpx2   ## For references, Davies+18's version used 18 and 4
    minpix = 10

    
    ### The actual fitting
    C = np.zeros(3*nfit)   ## Reserving space; this size is a requirement for the continuum fitter

    fit_sdss_cont(wave, flux, sigma, nfit, z_in, deltapix1, deltapix2, minpix, slopethresh, fluxthresh, C, Fscale) # The fit
    

    wave_array = wave_rest_array*(1+z_in)
    contt = get_cont(C,wave_array) ## This is how the resulting fit is extracted - C holds the cubic parameters resulting 
                             ## from the fitting, and get_cont interpolates that cubic onto any wavelength array passed as an input.
                             ## Now contt contains the smoothed continuum.

    fval = np.interp(wave_norm,wave_rest_array,contt)  # Normalising output, required

    if fval==0:             # check for fval=0
        print('interpolated value is 0; flux values and fit set to 0')
        flux_norm = np.zeros(len(wave_rest))        
        sigma_norm = np.zeros(len(wave_rest))       
        contblue = np.zeros(len(wave_array))        
    else:
        flux_norm = flux/fval           # normalised flux
        sigma_norm = sigma/fval         # normalised errors
        contblue=contt/fval             # normalised spline fit


    if np.isnan(contblue).any():        # Check for nans
        print('nan in fit, fit replaced with zeros')
        contblue = np.zeros(len(wave_array))


    spec_fit = [wave_rest_array,contblue]       # Spectrum fit in rest-frame
    spec_norm = [wave_rest,flux_norm]           # Normalised spectrum in rest-frame

    return spec_fit, spec_norm, sigma_norm


