import numpy as np
import glob
import os
import matplotlib.pyplot as plt

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

dpix= 25
ftype= 'dpx{}'.format(dpix)
ntype= 'norm'


xpath = '/media/bartosz/USB STICK/highz_data/'
xmeta = np.loadtxt(xpath+'meta_data_v2.txt', dtype='str')
xnames = xmeta[:,0]

xZ = xmeta[:,1]
#wl1r = [float(wl) for wl in xmeta[:,2]]

w1l = 13450     # water absorption observed WL
w1r = 14250


Wave = np.arange(1000,3000,0.5)               # wavelength array of fits is the same for all

Fits, Idxs = [], []
qso_list=[]
for h,sf in enumerate(spec_folder):
    spath = dpath + sf
    fpath = spath + 'fits/'                        # path to fits
    npath = spath + 'normed/'                      # path to normalised spectra

    flist = natsorted([f for f in glob.glob(fpath + '*.txt')])
    qso_list = np.append(qso_list,flist)
    meta_ext = np.loadtxt(dpath + meta_data_ext[h] + '.txt', dtype='str')
    idx_meta = [int(m) for m in meta_ext[:,0]]
    z_meta = [float(z) for z in meta_ext[:,1]]


    for i,s in enumerate(flist):

        # extract file name
        sfile_name = Path(s).stem                
        sfile_info = sfile_name.split('_')
        
        ffile = sfile_name + '.txt'
        nfile = sfile_name + '_' + ntype + '.txt'

        idx = int(sfile_info[0])


        fqso = np.loadtxt(fpath + ffile)
        fit = fqso[:,1]    

    
        ## POWER LAW
        alpha = alphas[pl_idxs.index(idx)]
        beta = betas[pl_idxs.index(idx)]

        # SUBTRACTION
        fit = fit - power_law(Wave,alpha,beta)
#
#        # DIVISION 
#        fit = fit / power_law(Wave,alpha,beta)


        # prepare fits in array for NN search
        Idxs.append(idx)
        Fits.append(fit)

        print(i)


#for j,xqso in enumerate(xnames):
for j in np.array([37,38,39]):
    #j = 9
    xqso = xnames[j]
    z = float(xZ[j])

    chosen_wave = Wave
    wl_a = 1250
    wl_b = w1l/(1+z)
    wl_c = w1r/(1+z)
    wl_d = 2250
    aint = chosen_wave > wl_a
    bint = chosen_wave < wl_b
    cint = np.logical_and(aint,bint)
    dint = np.logical_and(chosen_wave>wl_c, chosen_wave<wl_d)
    nint = np.logical_or(cint,dint)


    X = [row[nint] for row in Fits]

    nbrs = NearestNeighbors(n_neighbors=6, algorithm='kd_tree', metric='euclidean', p=2).fit(X)

    #NN_meta_idx = []
    #NN_d = []
    result=[]
    for i,qso in enumerate(qso_list):
        chosen_fit = Fits[i]    

        #set wavelength interval for nearest neighbor search

        distances, indices = nbrs.kneighbors([chosen_fit[nint]])

        nbrs_idx = indices[0]
        nbrs_d = distances[0]

        d_mean = np.mean(nbrs_d[1:])  

        #result.append([Idxs[i],*meta_idx,*nbrs_d])
        result.append([Idxs[i],d_mean])
        print(xqso,i)


    with open(xpath + 'eBOSS_NN_trim/pl_corr/{}_NN.txt'.format(xqso),'w') as nnsave:
        for x in result:
            #nnsave.write('{} \t {} \t{} \t{} \t{} \t{} \t{} \t {} \t{} \t{} \t{} \t{} \t{} \n'.format(*x))
            nnsave.write('{} \t {} \n'.format(*x))
    nnsave.close()
    print(xqso)

