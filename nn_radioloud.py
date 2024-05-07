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

w1l = 13450     # water absorption observed WL
w1r = 14250
w2l = 18200
w2r = 19400

# Load eBOSS fits and store them in array
Fits, Wave, Idxs = [], [], []
excl=0
for i,s in enumerate(slist):

    sfile_name = Path(s).stem                
    sfile_info = sfile_name.split('_')
    
    ffile = sfile_name + '.txt'
    nfile = sfile_name + '_' + ntype + '.txt'

    
    fqso = np.loadtxt(fpath + ffile)
    wave_rest_array = fqso[:,0]
    fit = fqso[:,1]

    if fit[1600]==0:
        excl+=1
        continue

    fit = fit/fit[1600]    


    Fits.append(fit)
    Wave.append(wave_rest_array)

    idx = sfile_info[0]
    nspec = sfile_info[2]
    Idxs.append(idx)
    print(i)

print(excl)

# high-z data meta
xpath = '/media/bartosz/USB STICK/added_sample_XIE/'
xmeta = np.loadtxt(xpath+'notes', dtype='str')
xnames = xmeta[:,0]
xZ = xmeta[:,1]

xspath = xpath + 'fits/'
xlist = [glob.glob(xspath + '{}_*.txt'.format(name))[0] for name in xnames]
print(xlist)

starts = [1379.0, 1496.0, 1568.0, 1439.0, 1512.5, 1550.5, 1496.5, 1524.0, 1572.5]
w1l = [1970,2125,2200,2025,2125,2150,2131.5,2150,2225]
w1r = [14885.019969102386, 14885.570505437146, 14884.977639106755, 14887.090558996453, 14887.879568414128, 14892.260472995786, 14905.148993983797, 14912.29986657092, 14894.839438404864]
result=[]
for i,f in enumerate(xlist[:-1]):
    if i == 0:
        continue
    z_in = float(xZ[i])
    
    chosen_wave = Wave[0]
    wl_a = starts[i]    #1250
    wl_b = w1l[i] #13450 / (1+z_in)
    wl_c = w1r[i] / (1+z_in) #14875.355378717572       #14250 / (1+z_in)
    wl_d = 18200 / (1+z_in)
    wl_e = 19400 / (1+z_in)
    wl_f = 3000 #2828   #3000
    abint = np.logical_and(chosen_wave > wl_a,chosen_wave < wl_b)
    cdint = np.logical_and(chosen_wave > wl_c,chosen_wave < wl_d)
    efint = np.logical_and(chosen_wave > wl_e,chosen_wave < wl_f)
    adint = np.logical_or(abint,cdint)
    #cfint = np.logical_and(chosen_wave > wl_c,chosen_wave < wl_f)
    nint = np.logical_or(adint,efint) #np.logical_or(abint,cfint) 
    
    X = [row[nint] for row in Fits]

    nbrs = NearestNeighbors(n_neighbors=50, algorithm='kd_tree', metric='euclidean', p=2).fit(X)

    qname = Path(f).stem.split('_')[0]

    fit = np.loadtxt(f)

    wl = fit[:,0]
    flux = fit[:,1]

    distances, indices = nbrs.kneighbors([flux[nint]])      #find NNs of fit

    nbrs_idx = indices[0]
    nbrs_d = np.array([str(d) for d in distances[0]])

    meta_idx = [str(Idxs[j]) for j in nbrs_idx]      #idx of NN in eBOSS meta table

    result.append([str(qname),*meta_idx,*nbrs_d])
    print(f)

with open(xpath + 'rl_NN_full_sc_v2.txt','w') as nnsave:
    for x in result:
        nnsave.write('\t'.join(x) + '\n')
        #nnsave.write('{} \t {} \t{} \t{} \t{} \t{} \t {} \t{} \t{} \t{} \t{} \n'.format(*x))
nnsave.close()
