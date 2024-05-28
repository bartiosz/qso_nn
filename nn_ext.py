import numpy as np
import glob

from sklearn.neighbors import NearestNeighbors

from pathlib import Path

from natsort import natsorted

import itertools

def power_law(x,alpha,beta):
    return alpha*x**beta

dpath= '/media/bartosz/USB STICK/BOSS_DR14_ext/'      # path to data
spec_folder = ['spectra_07_2', 'spectra_2_3', 'spectra_3_4']


# power law fit parameters
power_law = np.loadtxt(dpath + 'power_law_fits.txt', dtype='str')
pl_idxs = np.array([int(i) for i in power_law[:,0]])
alphas = np.array([float(a) for a in power_law[:,2]])
betas = np.array([float(b) for b in power_law[:,4]])

def power_law(x,alpha,beta):
    return alpha*x**beta


# select fit resolution
dpix= 25
ftype= 'dpx{}'.format(dpix)
ntype= 'norm'

Wave = np.arange(1000,3000,0.5)               # wavelength array of fits is the same for all
Fits, Idxs = [], []

for sf in spec_folder:
    spath = dpath + sf
    fpath = spath + 'fits/'                        # path to fits
    npath = spath + 'normed/'                      # path to normalised spectra

    flist = natsorted([f for f in glob.glob(fpath + '*.txt')])


    for i,s in enumerate(flist):

        # extract file name
        sfile_name = Path(s).stem                
        sfile_info = sfile_name.split('_')
        
        ffile = sfile_name + '.txt'
        nfile = sfile_name + '_' + ntype + '.txt'

        fqso = np.loadtxt(fpath + ffile)
        fit = fqso[:,1]    

    
#        ## POWER LAW
#        alpha = alphas[pl_idxs.index(idx)]
#        beta = betas[pl_idxs.index(idx)]
#
#        # SUBTRACTION
#        fit = fit - power_law(Wave,alpha,beta)
#
#        # DIVISION 
#        fit = fit / power_law(Wave,alpha,beta)


        # prepare fits in array for NN search
        Fits.append(fit)

        idx = sfile_info[0]
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

# power law fits of highz quasars
pl = np.loadtxt(xpath + 'power_law_fits.txt', dtype='str')
alphas = np.array([float(a) for a in pl[:,1])
betas = np.array([float(b) for b in pl[:,3])

# find best neighbors for every quasar in list
result=[]
for i,qso in enumerate(xlist):

    chosen_wave = Wave
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

    # ## POWER LAW
    # alpha = alphas[pl[:,0].index(qso)]
    # beta = betas[pl[:,0].index(qso)]
    
    # # SUBTRACTION
    # flux = flux + power_law(wl,alpha,beta)

    # # Division
    # flux = flux/power_law(wl,alpha,beta)
    
    # calculate NNs and their distance
    distances, indices = nbrs.kneighbors([flux[nint]])

    nbrs_idx = indices[0]
    nbrs_d = distances[0]
    
    meta_idx = [Idxs[j] for j in nbrs_idx]    

    result.append([qname,*meta_idx,*nbrs_d])

    print(i)


# save results
with open(xpath + 'highZ_NN_full_ext.txt','w') as nnsave:
    for x in result:
        nnsave.write('{} \t {} \t{} \t{} \t{} \t{} \t {} \t{} \t{} \t{} \t{} \n'.format(*x))
nnsave.close()







