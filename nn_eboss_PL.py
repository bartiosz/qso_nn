import numpy as np
import glob

from sklearn.neighbors import NearestNeighbors

from pathlib import Path

from natsort import natsorted

import itertools

dpath= '/media/bartosz/USB STICK/BOSS_DR14_ext/'      # path to data
spec_folder = ['spectra_07_2/', 'spectra_2_3/', 'spectra_3_4/']
meta_data_ext = ['meta_data_ext_07_2','meta_data_ext_2_3','meta_data_ext_3_4']

# power law fit parameters
pl = np.loadtxt(dpath + 'power_law_fits_spline_rsq_cl_bound.txt', dtype='str')
pl_idxs = [int(i) for i in pl[:,0]]
alphas = [float(a) for a in pl[:,2]]
betas = [float(b) for b in pl[:,4]]

def power_law(x,alpha,beta):
    return alpha*x**beta


# select fit resolution
dpix= 25
ftype= 'dpx{}'.format(dpix)
ntype= 'norm'

Wave = np.arange(1000,3000,0.5)               # wavelength array of fits is the same for all
Fits, Idxs = [], []

for h,sf in enumerate(spec_folder):
    spath = dpath + sf
    fpath = spath + 'fits/'                        # path to fits
    npath = spath + 'normed/'                      # path to normalised spectra

    flist = natsorted([f for f in glob.glob(fpath + '*.txt')])
    meta_ext = np.loadtxt(dpath + meta_data_ext[h] + '.txt', dtype='str')
    idx_meta = [int(m) for m in meta_ext[:,0]]
    z_meta = [float(z) for z in meta_ext[:,1]]


    for i,s in enumerate(flist):

        # extract file name
        sfile_name = Path(s).stem                
        sfile_info = sfile_name.split('_')
        
        ffile = sfile_name + '.txt'
        nfile = sfile_name + '_' + ntype + '.txt'

        fqso = np.loadtxt(fpath + ffile)
        fit = fqso[:,1]    

        idx = int(sfile_info[0])

        ## POWER LAW
        alpha = alphas[pl_idxs.index(idx)]
        beta = betas[pl_idxs.index(idx)]

        # SUBTRACTION
        fit = fit - power_law(Wave,alpha,beta)

        # prepare fits in array for NN search
        Fits.append(fit)
        Idxs.append(idx)
        print(i)


# define analysis WL range
chosen_wave = Wave
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
result=[]
for h,sf in enumerate(spec_folder):
    spath = dpath + sf
    fpath = spath + 'fits/'                        # path to fits
    npath = spath + 'normed/'                      # path to normalised spectra

    qso_list = natsorted([f for f in glob.glob(fpath + '*.txt')])
    
    for i,qso in enumerate(qso_list):

        # determine file name of spectrum
        chosen_file = Path(qso).stem
        chosen_info = chosen_file.split('_')

        chosen = int(chosen_info[0])

        # load spectrum and fit
        chosen_idx = Idxs.index(chosen)
        chosen_fit = Fits[chosen_idx]
        chosen_wave = Wave

        # calculate NNs and their distance
        distances, indices = nbrs.kneighbors([chosen_fit[nint]])

        nbrs_idx = indices[0]
        nbrs_d = distances[0]
        
        meta_idx = [Idxs[j] for j in nbrs_idx]    

        result.append([Idxs[i],*meta_idx,*nbrs_d])

        print(i)


# save results
with open(dpath + 'nearest_neighbors_PLcorr.txt','w') as nnsave:
    for x in result:
        nnsave.write('{} \t \t {} \t{} \t{} \t{} \t{} \t{} \t\t {} \t{} \t{} \t{} \t{} \t{} \n'.format(*x))
nnsave.close()







