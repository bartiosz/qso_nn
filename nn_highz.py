import numpy as np
import glob
import os
import matplotlib.pyplot as plt

from sklearn.neighbors import NearestNeighbors

from pathlib import Path

from natsort import natsorted

import itertools

dpath= '/media/bartosz/USB STICK/BOSS_DR14/'
spath= dpath + 'spectra/'
fpath= dpath + 'fits/'
npath= dpath + 'normed/'

dpix= 25
ftype= 'dpx{}'.format(dpix)
ntype= 'norm'


slist = natsorted([s for s in glob.glob(fpath + '*.txt')])


# Load eBOSS fits and store them in array
Fits, Wave, Idxs = [], [], []
for i,s in enumerate(slist):

    sfile_name = Path(s).stem                
    sfile_info = sfile_name.split('_')
    
    ffile = sfile_name + '.txt'
    nfile = sfile_name + '_' + ntype + '.txt'

    
    fqso = np.loadtxt(fpath + ffile)
    wave_rest_array = fqso[:,0]
    fit = fqso[:,1]    


    Fits.append(fit)
    Wave.append(wave_rest_array)

    idx = sfile_info[0]
    nspec = sfile_info[2]
    Idxs.append(idx)
    print(i)


# high-z data meta
xpath = '/media/bartosz/USB STICK/highz_data/'
xmeta = np.loadtxt(xpath+'meta_data_v2.txt', dtype='str')
xnames = xmeta[:,0]
xZ = xmeta[:,1]
wl1r = [float(wl) for wl in xmeta[:,3]]
wl2r = [float(wl) for wl in xmeta[:,4]]


xspath = xpath + 'fits/'
xlist = [glob.glob(xspath + '{}_*.txt'.format(name))[0] for name in xnames]
print(xlist)


result=[]
for i,f in enumerate(xlist):

    chosen_wave = Wave[0]
    wl_a = 1250
    wl_b = wl1r[i]
    wl_c = wl2r[i]
    wl_d = 2250
    abint = np.logical_and(chosen_wave > wl_a,chosen_wave < wl_b)
    cdint = np.logical_and(chosen_wave > wl_c,chosen_wave < 2250)
    nint = np.logical_or(abint,cdint)
    
    X = [row[nint] for row in Fits]

    nbrs = NearestNeighbors(n_neighbors=5, algorithm='ball_tree', metric='euclidean', p=2).fit(X)

    qname = Path(f).stem.split('_')[0]

    fit = np.loadtxt(f)

    wl = fit[:,0]
    flux = fit[:,1]

    distances, indices = nbrs.kneighbors([flux[nint]])      #find NNs of fit

    nbrs_idx = indices[0]
    nbrs_d = distances[0]

    meta_idx = [Idxs[j] for j in nbrs_idx]      #idx of NN in eBOSS meta table

    result.append([qname,*meta_idx,*nbrs_d])
    print(f)


with open(xpath + 'highZ_NN_full_sc.txt','w') as nnsave:
    for x in result:
        nnsave.write('{} \t {} \t{} \t{} \t{} \t{} \t {} \t{} \t{} \t{} \t{} \n'.format(*x))
nnsave.close()
