import os
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord

from specdb import query_catalog as spqcat
from specdb import utils as spdbu
from specdb.specdb import SpecDB, IgmSpec
from specdb import specdb as sdbsdb

from pyigm.surveys.llssurvey import LLSSurvey
import matplotlib.pyplot as plt

import importlib
import random

import pickle

# select redshift range snd SNR threshold
z_min = 3
z_max = 4
SNR_min = 5


# results path
rpath= '/media/bartosz/USB STICK/BOSS_DR14_ext/'          

# data file    
db_file = '/media/bartosz/USB STICK/BOSS_DR14/IGMspec_DB_v03.1.hdf5'         

importlib.reload(sdbsdb)
sdb = sdbsdb.SpecDB(db_file=db_file)


# meta data
all_meta = sdb['BOSS_DR14'].meta            


# filter for redshifts
index_z = [i for i,entry in enumerate(all_meta) if z_min<=entry[6]<z_max]


# load and filter SNRs
snr_file = rpath + 'dr14_median_snr.pckl'
median_snr = pickle.load(open(snr_file,'rb'))

index_snr = np.array([i for i,snr in enumerate(median_snr) if snr>=SNR_min])
meta_snr = [median_snr[i] for i in index_snr]


# crossmatch redshift and SNR filtered lists
crossmatch = set(index_snr).intersection(index_z)
index_final = list(crossmatch)
index_final.sort()


##SNR and redshift of selected objects
snr_final = [median_snr[idx] for idx in index_final]
z_final = [all_meta[idx][6] for idx in index_final]


# meta data of selected objects
meta_final = sdb['BOSS_DR14'].meta_from_ids(np.array(index_final))


# write meta data file
meta_file = open(rpath + 'meta_data_ext.txt', 'w')

# load and save spectra of selected objects
e = 0               # e to count exceptions (not able to save)
for i,idx in enumerate(index_final):
    print(idx,  '(', i, ' of ', len(index_final), ')')

    # load quasar properties
    meta = meta_final[i]
    snr = snr_final[i]
    redshift = meta[6]
    plate,mjd,fiber = meta[3],meta[4],meta[5]
    specid = meta[-1]

    try:
        spec, meta = sdb.spectra_from_ID(idx)

    except Exception as error:
        saved = False
        err = type(error).__name__
        print('An exception occured:', err)
        e+=1
    
    else:        
        # select latest spectrum of object
        n = 0
        spec.select = n
        wl, flux, sig = spec.wavelength, spec.flux, spec.sig

        file_save='{}_spec-{}-{}-{}_{}'.format(idx,plate,mjd,fiber,n)
        with open(rpath + 'spectra/' + file_save + '.txt','w') as f:
            lst = [wl.value,flux,sig]
            for x in zip(*lst):
                f.write('{} \t {} \t {} \n'.format(*x))
        f.close()
    
        saved = True
        err = '-'

    # save meta
    meta_file.write('{} \t {} \t {} \t {} \t {} \t {} \n'.format(idx,redshift,snr,specid,saved,err))
    

meta_file.close()
print('number of exceptions: ', e)
