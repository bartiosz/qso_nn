from spline_fit_v2 import spline_fit
from funcs import get_z

import numpy as np

import glob

from specdb.specdb import SpecDB, IgmSpec
from specdb import specdb as sdbsdb

import importlib

from pathlib import Path

import itertools

# redshift range
#spec_folder = 'spectra__07/'
#spec_folder = 'spectra_07_2/'
#spec_folder = 'spectra_2_3/'
#spec_folder = 'spectra_3_4/'


db_file = '/media/bartosz/USB STICK/BOSS_DR14/IGMspec_DB_v03.1.hdf5'       # Path to catalog data
dpath = '/media/bartosz/USB STICK/BOSS_DR14_ext/' + spec_folder            # Path to spectra data
rpath = dpath + 'fits/'                                                    # where to save fit
npath = dpath + 'normed/'                                                  # where to save normalised spectra

# Load meta data for redshift
importlib.reload(sdbsdb)
sdb = sdbsdb.SpecDB(db_file=db_file)
all_meta = sdb['BOSS_DR14'].meta


# List of spectra to fit
speclist = [f for f in glob.glob(dpath + '*.txt')]


for i,f in enumerate(speclist):
    
    # Get quasar ID
    file_name = Path(f).stem
    file_info = file_name.split('_')

    
    # Load redshift
    qso_idx = int(file_info[0])
    z_in = get_z(all_meta,qso_idx)

#    if z_in < 1.8:   # spline has issues fitting quasars with z < 1.8
#        continue
#    else:
#        print(z_in)


    # Load spectrum data
    qso = np.loadtxt(f)
    wave = qso[:,0]
    flux = qso[:,1]
    sigma = qso[:,2]


    # Fit spline to spectrum
    deltapix2 = 25
    spec_fit, spec_norm, sigma_norm = spline_fit(wave,flux,sigma,z_in,deltapix2)


    # Save normalised spectrum
    norm_save='{}_norm'.format(file_name)
    with open(npath + norm_save + '.txt','w') as nsave:
        for x in itertools.zip_longest(*spec_norm,sigma_norm):
            nsave.write('{} \t {} \t {}\n'.format(*x))
    nsave.close()


    # Save normalised fit
    fit_save='{}_dpx{}'.format(file_name,deltapix2)
    with open(rpath + fit_save + '.txt','w') as fsave:
        for x in itertools.zip_longest(*spec_fit):
            fsave.write('{} \t {}\n'.format(*x))
    fsave.close()

    print(file_name, ', ', i, ' of ', len(speclist))





