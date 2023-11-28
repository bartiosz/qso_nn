import numpy as np
import glob
import os
import matplotlib.pyplot as plt

from sklearn.neighbors import NearestNeighbors

from pathlib import Path

from natsort import natsorted

import itertools

dpath= '/media/bartosz/Volume/BOSS_DR14/data/'
spath= dpath + 'spectra/'
fpath= dpath + 'fits/dpx25/'
npath= dpath + 'normed/dpx25/'

dpix= 25
ftype= 'dpx{}'.format(dpix)
ntype= 'norm'


slist = natsorted([s for s in glob.glob(spath + '*_0.txt')])


Fits, Wave, Idxs = [], [], []
for i,s in enumerate(slist):

    sfile_name = Path(s).stem                
    sfile_info = sfile_name.split('_')
    
    ffile = sfile_name + '_' + ftype + '.txt'
    nfile = sfile_name + '_' + ntype + '.txt'

    
    fqso = np.loadtxt(fpath + ffile)
    wave_rest_array = fqso[:,0]
    fit = fqso[:,1]    
    
    #nqso = np.loadtxt(npath + nfile)
    #wave_rest = nqso[:,0]
    #flux_norm = nqso[:,1]
    #sigma_norm = nqso[:,2] 

    Fits.append(fit)
    Wave.append(wave_rest_array)

    idx = sfile_info[0]
    nspec = sfile_info[2]
    Idxs.append(idx)
    print(i)


xpath = '/media/bartosz/Volume/XQR30/data/'
xmeta = np.loadtxt(xpath+'meta.txt', dtype='str')
xnames = xmeta[:,0]
xZ = xmeta[:,1]
wl1r = [float(wl) for wl in xmeta[:,2]]


xspath = xpath + 'fits/'
xlist = [glob.glob(xspath + '{}_*.txt'.format(name))[0] for name in xnames]
print(xlist)
#xlist = [s for s in glob.glob(xspath + '*_dpx34.txt')]

result=[]
for i,f in enumerate(xlist):

    chosen_wave = Wave[0]
    wl_a = 1250
    wl_b = wl1r[i]
    aint = chosen_wave > wl_a
    bint = chosen_wave < wl_b
    nint = np.logical_and(aint,bint)


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


with open(xpath + 'highZ_NN_trim.txt','w') as nnsave:
    for x in result:
        nnsave.write('{} \t \t {} \t{} \t{} \t{} \t{} \t\t {} \t{} \t{} \t{} \t{} \n'.format(*x))
nnsave.close()
