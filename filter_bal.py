import numpy as np
import importlib
import glob
from funcs import plot_fit, get_meta, get_z, get_BI

import shutil

from specdb.specdb import SpecDB, IgmSpec
from pathlib import Path

from natsort import natsorted

dpath= '/media/bartosz/Volume/BOSS_DR14/data/'
spath= dpath + 'spectra/'
fpath= dpath + 'fits/dpx25/'
npath= dpath + 'normed/dpx25/'

igmsp = SpecDB('/media/bartosz/Volume/igmspec_data/DB/IGMspec_DB_v03.1.hdf5')
meta = igmsp['BOSS_DR14'].meta

slist = natsorted([s for s in glob.glob(spath + '*_0.txt')])
flist = natsorted([s for s in glob.glob(fpath + '*_0_dpx25.txt')])
nlist = natsorted([s for s in glob.glob(npath + '*_0_norm.txt')])

for i, s in enumerate(slist):
    sfile_name = Path(s).stem                
    sfile_info = sfile_name.split('_')

    idx = int(sfile_info[0])

    BI = get_BI(meta,idx)

    if BI > 0:
        shutil.move(s, spath+'BALs/')
        shutil.move(flist[i], fpath+'BALs/')
        shutil.move(nlist[i], npath+'BALs/')

    print(i)
