from spline_fit_v2 import spline_fit

import numpy as np

import glob

import importlib

import os
from pathlib import Path

import itertools


dpath = '/media/bartosz/Volume/highz_data/'
spath = dpath + 'spectra/'
rpath = dpath + 'fits/'
npath = dpath + 'normed/'

# Spline sensitivity
deltapix2 = round(25 * 6.902467685076381)

# Load meta for redshift
meta = np.loadtxt(dpath + 'meta_data_v2.txt',dtype='str')
xnames = meta[:,0]
Z = meta[:,1]
wl1r = [float(wl) for wl in meta[:,3]]      #short rest-frame wavelength of water aborption range
wl2r = [float(wl) for wl in meta[:,4]]      #long ---''---


# List of spectra to fit
xlist = [glob.glob(spath + '{}.txt'.format(name))[0] for name in xnames]

for i,f in enumerate(xlist):
    
    # Get quasar ID
    file_name = Path(f).stem
    file_info = file_name.split('_')
    
    # Load redshift
    qso_idx = np.nonzero(xnames == file_info[0])[0][0]
    z_in = float(Z[qso_idx])

    
    # Load quasar spectrum
    qso = np.loadtxt(f)
    wave = qso[:,0]
    flux = qso[:,1]
    sigma = qso[:,2]


    # Fit spline to spectrum (include water absorption)
    spec_fit, spec_norm, sigma_norm = spline_fit(wave,flux,sigma,z_in,deltapix2)


    # mask water absorption
    wl1 = 13450
    wl2 = 14250
    water_mask_1 = wave < wl1
    water_mask_2 = wave > wl2
    water_mask = np.logical_or(water_mask_1,water_mask_2)
    print(water_mask)
    wave, flux, sigma = wave[water_mask], flux[water_mask], sigma[water_mask]


    # Fit spline to spectrum (exclude water absorption)
    spec_fit, spec_norm2, sigma_norm2 = spline_fit(wave,flux,sigma,z_in,deltapix2)


    # Save normalised spectrum
    norm_save='{}_full_norm'.format(file_name)
    with open(npath + norm_save + '.txt','w') as nsave:
        for x in itertools.zip_longest(*spec_norm,sigma_norm):
            nsave.write('{} \t {} \t {}\n'.format(*x))
    nsave.close()


    # save spectrum fit
    fit_save='{}_full_dpx{}'.format(file_name,deltapix2)
    with open(rpath + fit_save + '.txt','w') as fsave:
        for x in itertools.zip_longest(*spec_fit):
            fsave.write('{} \t {}\n'.format(*x))
    fsave.close()


    print(file_name)





