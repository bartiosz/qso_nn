import numpy as np
import matplotlib.pyplot as plt
from fit_cont._cont_src import fit_sdss_cont

import glob

from sklearn.neighbors import NearestNeighbors

from specdb import query_catalog as spqcat
#from specdb import interface_db as spgidb
from specdb import utils as spdbu
from specdb.specdb import SpecDB, IgmSpec
from specdb import specdb as sdbsdb

import importlib

import os
from pathlib import Path

import itertools


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

dpath= '/media/bartosz/Volume/BOSS_DR14/data/spectra/'
db_file = '/media/bartosz/Volume/igmspec_data/DB/IGMspec_DB_v03.1.hdf5'
rpath = '/media/bartosz/Volume/BOSS_DR14/data/fits/'
npath = '/media/bartosz/Volume/BOSS_DR14/data/normed/'

#load meta data for redshift
importlib.reload(sdbsdb)
sdb = sdbsdb.SpecDB(db_file=db_file)
all_meta = sdb['BOSS_DR14'].meta


# Automated continuum fitter parameters. Built-in params for steepness and overall flux adjustment.
# Avoid changing these.
slopethresh = 0.033 
fluxthresh = 0.99
Fscale = 1.0

# Normalisation wavelength
wave_norm = 1290.0


#z = [2.865,2.656,2.254,2.555,2.988,2.324]
fits=[]
fvals=[]
mylist = [f for f in glob.glob(dpath + '*_0.txt')]


wave_rest_array = np.arange(1000,3000,0.5)

#interval to match neighbors
a_int = 500
b_int = 2500

a_wl = wave_rest_array[a_int]
b_wl = wave_rest_array[b_int]

#fig, ax = plt.subplots(figsize=(28,7))
for i,f in enumerate(mylist):
    #if i==100: break;
    
    file_name = Path(f).stem
    file_info = file_name.split('_')
    
    #load redshift
    qso_idx = int(file_info[0])
    qso_meta = all_meta[qso_idx]
    z_in = qso_meta[6]

    qso = np.loadtxt(f)
    wave = qso[:,0]
    flux = qso[:,1]
    sigma = qso[:,2]
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
    deltapix2 = 25   ## For references, Davies+18's version used 18 and 4
    minpix = 10

    
    ### The actual fitting
    C = np.zeros(3*nfit)   ## Reserving space; this size is a requirement for the continuum fitter

    fit_sdss_cont(wave, flux, sigma, nfit, z_in, deltapix1, deltapix2, minpix, slopethresh, fluxthresh, C, Fscale) # The fit
    

    wave_array = wave_rest_array*(1+z_in)
    contt = get_cont(C,wave_array) ## This is how the resulting fit is extracted - C holds the cubic parameters resulting 
                             ## from the fitting, and get_cont interpolates that cubic onto any wavelength array passed as an input.
                             ## Now contt contains the smoothed continuum.

    fval = np.interp(wave_norm,wave_rest_array,contt)  # Normalising output, required
    fvals.append(fval)
    if fval==0:
        print(file_name, 'interpolated value is 0; flux values and fit set to 0')
        flux_norm = np.zeros(len(wave_rest))
        sigma_norm = np.zeros(len(wave_rest))
        contblue = np.zeros(len(wave_array))
    else:
        flux_norm = flux/fval
        sigma_norm = sigma/fval
        contblue=contt/fval



    if np.isnan(contblue).any():
        print(file_name, 'nan in fit, fit replaced with zeros')
        contblue = np.zeros(len(wave_array))


    norm_save='{}_norm'.format(file_name)
    with open(npath + norm_save + '.txt','w') as nsave:
        for x in itertools.zip_longest(wave_rest,flux_norm,sigma_norm):
            nsave.write('{} \t {} \t {}\n'.format(*x))
    nsave.close()


    fit_save='{}_dpx{}'.format(file_name,deltapix2)
    with open(rpath + fit_save + '.txt','w') as fsave:
        for x in itertools.zip_longest(wave_rest_array,contblue):
            fsave.write('{} \t {}\n'.format(*x))
    fsave.close()

    print(i)

#fval_save = 'fvals.txt'
#s = open(rpath + fval_save, 'w')
#s.write('{}'.format(fvals))
#s.close()




