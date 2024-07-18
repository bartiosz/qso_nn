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

        idx = int(sfile_info[0])


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
        Idxs.append(idx)
        Fits.append(fit)

        print(i)


# high-z data meta
xpath = '/media/bartosz/USB STICK/highz_data/'
xmeta = np.loadtxt(xpath+'meta_data_v2.txt', dtype='str')
xnames = xmeta[:,0]
xZ = xmeta[:,1]
wl1r = [float(wl) for wl in xmeta[:,3]]         # water absorption range
wl2r = [float(wl) for wl in xmeta[:,4]]

# quasar list
xspath = xpath + 'fits/'
xlist = [glob.glob(xspath + '{}_*.txt'.format(name))[0] for name in xnames]
print(xlist)

# load power law fit for high z quasars
pl = np.loadtxt('/media/bartosz/USB STICK/highz_data/power_law_fits_spline_rsq_cl_bound.txt',dtype='str')
pl_name = pl[:,0]
alphas = [float(a) for a in pl[:,1]]
betas = [float(b) for b in pl[:,3]]


result=[]
for i,f in enumerate(xlist):

    # define analysis WL range
    chosen_wave = Wave
    wl_a = 1250
    wl_b = wl1r[i]
    wl_c = wl2r[i]
    wl_d = 2250
    abint = np.logical_and(chosen_wave > wl_a,chosen_wave < wl_b)
    cdint = np.logical_and(chosen_wave > wl_c,chosen_wave < wl_d)
    nint = np.logical_or(abint,cdint)
    
    # mask for analysis range
    X = [row[nint] for row in Fits]
    
    # ball tree arrangement of data
    nbrs = NearestNeighbors(n_neighbors=5, algorithm='ball_tree', metric='euclidean', p=2).fit(X)

    qname = Path(f).stem.split('_')[0]

    # load spectrum and fit
    fit = np.loadtxt(f)

    wl = fit[:,0]
    flux = fit[:,1]

#
#    # power law 
#    pl_idx = np.where(pl_name == qname)[0][0]
#    alpha = alphas[pl_idx]
#    beta = betas[pl_idx]
#
#    # subtraction
#    flux = flux - power_law(Wave,alpha,beta)
#
#    # division
#    flux = flux / power_law(Wave,alpha,beta)
    

    # calculate NNs and their distance
    distances, indices = nbrs.kneighbors([flux[nint]])      #find NNs of fit

    nbrs_idx = indices[0]
    nbrs_d = distances[0]

    meta_idx = [Idxs[j] for j in nbrs_idx]      #idx of NN in eBOSS meta table

    result.append([qname,*meta_idx,*nbrs_d])
    print(f)


with open(xpath + 'highZ_NN_uncorr.txt','w') as nnsave:
    for x in result:
        nnsave.write('{} \t {} \t{} \t{} \t{} \t{} \t {} \t{} \t{} \t{} \t{} \n'.format(*x))
nnsave.close()




