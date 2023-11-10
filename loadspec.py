import os
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord

from specdb import query_catalog as spqcat
#from specdb import interface_db as spgidb
from specdb import utils as spdbu
from specdb.specdb import SpecDB, IgmSpec
from specdb import specdb as sdbsdb

from pyigm.surveys.llssurvey import LLSSurvey
import matplotlib.pyplot as plt

import importlib
import random

import pickle

rpath= '/media/bartosz/Volume/BOSS_DR14/'
db_file = '/media/bartosz/Volume/igmspec_data/DB/IGMspec_DB_v03.1.hdf5'

importlib.reload(sdbsdb)
sdb = sdbsdb.SpecDB(db_file=db_file)

all_meta = sdb['BOSS_DR14'].meta

index_z = [i for i,entry in enumerate(all_meta) if 2<=entry[6]<=3]


snr_file = rpath + 'dr14_median_snr.pckl'
median_snr = pickle.load(open(snr_file,'rb'))

index_snr = np.array([i for i,snr in enumerate(median_snr) if snr>=20])
meta_snr = [median_snr[i] for i in index_snr]

crossmatch = set(index_snr).intersection(index_z)
index_final = list(crossmatch)
index_final.sort()

snr_final = [median_snr[idx] for idx in index_final]
z_final = [all_meta[idx][6] for idx in index_final]


meta_final = sdb['BOSS_DR14'].meta_from_ids(np.array(index_final))


for i,idx in enumerate(index_final):
    spec, meta = sdb.spectra_from_ID(idx)

    meta = meta_final[i]
    plate,mjd,fiber=meta[3],meta[4],meta[5]

    n = 0
    spec.select = n
    wl, flux, sig = spec.wavelength, spec.flux, spec.sig

    file_save='{}_spec-{}-{}-{}_{}'.format(idx,plate,mjd,fiber,n)
    with open(rpath + 'data/' + file_save + '.txt','w') as f:
        lst = [wl.value,flux,sig]
        for x in zip(*lst):
            f.write('{} \t {} \t {} \n'.format(*x))
    f.close()
