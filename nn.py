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


chosen_wave = Wave[0]
wl_a = 1250
wl_b = 2250
aint = chosen_wave > wl_a
bint = chosen_wave < wl_b
nint = np.logical_and(aint,bint)


X = [row[nint] for row in Fits]

nbrs = NearestNeighbors(n_neighbors=6, algorithm='ball_tree', metric='euclidean', p=2).fit(X)

qso_list = slist

#NN_meta_idx = []
#NN_d = []
result=[]
for i,qso in enumerate(qso_list):

    chosen_file = Path(qso).stem
    chosen_info = chosen_file.split('_')

    chosen = chosen_info[0]

    chosen_idx = Idxs.index(chosen)
    chosen_fit = Fits[chosen_idx]
    chosen_wave = Wave[chosen_idx]

    #sfile_name = Path(slist[chosen_idx]).stem                

    #nfile = sfile_name + '_' + ntype + '.txt'
    #chosen_qso = np.loadtxt(npath + nfile)
    #wave_rest = chosen_qso[:,0]
    #flux_norm = chosen_qso[:,1]
    #sigma_norm = chosen_qso[:,2] 


    #set wavelength interval for nearest neighbor search

    distances, indices = nbrs.kneighbors([chosen_fit[nint]])

    nbrs_idx = indices[0]
    nbrs_d = distances[0]
    
    meta_idx = [Idxs[j] for j in nbrs_idx]    

    result.append([Idxs[i],*meta_idx,*nbrs_d])

    #NN_meta_idx.append(meta_idx)
    #NN_d.append(nbrs_d)

    print(i)



with open(dpath + 'nearest_neighbors_nonBALs_2.txt','w') as nnsave:
    for x in result:
        nnsave.write('{} \t \t {} \t{} \t{} \t{} \t{} \t{} \t\t {} \t{} \t{} \t{} \t{} \t{} \n'.format(*x))
nnsave.close()







