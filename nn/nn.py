import numpy as np
import glob

from sklearn.neighbors import NearestNeighbors

from pathlib import Path

from natsort import natsorted

import itertools

dpath= '/media/bartosz/Volume/BOSS_DR14/'      # path to data
spath= dpath + 'spectra/'                           # path to spectra  
fpath= dpath + 'fits/'                        # path to fits
npath= dpath + 'normed/'                      # path to normalised spectra


# select fit resolution
dpix= 25
ftype= 'dpx{}'.format(dpix)
ntype= 'norm'


slist = natsorted([s for s in glob.glob(spath + '*_0.txt')])


Fits, Wave, Idxs = [], [], []
for i,s in enumerate(slist):

    # extract file name
    sfile_name = Path(s).stem                
    sfile_info = sfile_name.split('_')
    
    ffile = sfile_name + '_' + ftype + '.txt'
    nfile = sfile_name + '_' + ntype + '.txt'

    fqso = np.loadtxt(fpath + ffile)
    wave_rest_array = fqso[:,0]
    fit = fqso[:,1]    

    # prepare fits in array for NN search
    Fits.append(fit)
    Wave.append(wave_rest_array)

    idx = sfile_info[0]
    Idxs.append(idx)
    print(i)


# define analysis WL range
chosen_wave = Wave[0]
wl_a = 1250
wl_b = 2250
aint = chosen_wave > wl_a
bint = chosen_wave < wl_b
nint = np.logical_and(aint,bint)


# mask for analysis range
X = [row[nint] for row in Fits]


# ball tree arrangement of data
nbrs = NearestNeighbors(n_neighbors=6, algorithm='ball_tree', metric='euclidean', p=2).fit(X)


# find best neighbors for every quasar in list
qso_list = slist
result=[]
for i,qso in enumerate(qso_list):

    # determine file name of spectrum
    chosen_file = Path(qso).stem
    chosen_info = chosen_file.split('_')

    chosen = chosen_info[0]

    # load spectrum and fit
    chosen_idx = Idxs.index(chosen)
    chosen_fit = Fits[chosen_idx]
    chosen_wave = Wave[chosen_idx]

    # calculate NNs and their distance
    distances, indices = nbrs.kneighbors([chosen_fit[nint]])

    nbrs_idx = indices[0]
    nbrs_d = distances[0]
    
    meta_idx = [Idxs[j] for j in nbrs_idx]    

    result.append([Idxs[i],*meta_idx,*nbrs_d])

    print(i)


# save results
with open(dpath + 'nearest_neighbors_snr10.txt','w') as nnsave:
    for x in result:
        nnsave.write('{} \t \t {} \t{} \t{} \t{} \t{} \t{} \t\t {} \t{} \t{} \t{} \t{} \t{} \n'.format(*x))
nnsave.close()







