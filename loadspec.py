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

def degtodegms(c):
    deg,cc=divmod(c,1)
    m,s=divmod(cc*60,1)
    return '{:.0f}° {:.0f}\' {:.1f}\'\''.format(deg,m,s)

def degtohms(c):
    h,cc=divmod(c/360*24,1)
    m,s=divmod(cc*60,1)
    return '{:.0f}h {:.0f}\' {:.1f}\'\''.format(h,m,s)

rpath= '/media/bartosz/Volume/BOSS_DR14/'
#rpath='/home/bartosz/Projects/BOSS/'
db_file = '/media/bartosz/Volume/igmspec_data/DB/IGMspec_DB_v03.1.hdf5'

importlib.reload(sdbsdb)
#igmsp = sdbsdb.IgmSpec(db_file='/media/bartosz/Volume/igmspec_data/DB/IGMspec_DB_v03.1.hdf5')
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
    #coords = SkyCoord(obj[0], obj[1],unit='deg')
    #meta_sing = sdb.meta_from_position(coords, 1*u.arcsec)
    #spec = sdb.spectra_from_meta(meta_sing)
    
    #for n in range(spec.nspec):
    n = 0
    spec.select = n
    wl, flux, sig = spec.wavelength, spec.flux, spec.sig

    file_save='{}_spec-{}-{}-{}_{}'.format(idx,plate,mjd,fiber,n)
    with open(rpath + 'data/' + file_save + '.txt','w') as f:
        lst = [wl.value,flux,sig]
        for x in zip(*lst):
            f.write('{} \t {} \t {} \n'.format(*x))
    f.close()

quit()
#######################





for i, obj in enumerate(meta):
    spec.select=i
    array_wl,array_flux,array_sig=spec.wavelength,spec.flux,spec.sig
    plate,mjd,fiber = obj[3],obj[4],obj[5]
    ra,dec,z=obj[0],obj[1],obj[6]
    file_save='spec-{}-{}-{}'.format(plate,mjd,fiber)
    f=open(rpath + 'data/' + file_save + '.txt','w')
    for j in range(len(array_wl)):
        wl,flux,sig=array_wl[j],array_flux[j],array_sig[j]
        f.write('{}\t{}\t{}\n'.format(wl,flux,sig))
    f.close()
    
    fig, ax = plt.subplots()

    ax.plot(array_wl,array_flux, lw=1,c='black')
    ax.set_title('{}'.format(file_save), fontsize=14)
    ax.set_xlabel('wavelength $[\AA]$',fontsize=14)
    ax.set_ylabel('flux',fontsize=14)
    fig.set_figwidth(15)
    info1 = '\n'.join((
        r'$RA= {}$'.format(degtohms(ra)),
        r'$DEC={}$'.format(degtodegms(dec)),
        r'$z={0:.3f}$'.format(z)))
    info2 = '\n'.join((
        r'plate: {}'.format(plate),
        r'mjd: {}'.format(mjd),
        r'fiber: {}'.format(fiber)))

    props = dict(boxstyle='square', facecolor='lightblue', alpha=0.5)
    ax.text(0.84, 0.95, info1, transform=ax.transAxes, fontsize=14,
            verticalalignment='top', bbox=props)
    ax.text(0.7, 0.95, info2, transform=ax.transAxes, fontsize=14,
            verticalalignment='top', bbox=props)
    plt.savefig(rpath + 'plots/{}.png'.format(file_save), format='png')
    plt.close()
    
    break
#spec.plot()
quit()

