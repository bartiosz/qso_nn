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

#qso_list = np.array(['12891_spec-6184-56267-236_0.txt','100632_spec-4501-55590-712_0.txt','103033_spec-4465-55858-962_0.txt','177263_spec-4637-55616-112_1.txt','199036_spec-5350-56009-475_1.txt'])
#qso_list = np.random.choice(np.array(slist),20)
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



with open(dpath + 'nearest_neighbors.txt','w') as nnsave:
    for x in result:
        nnsave.write('{} \t \t {} \t{} \t{} \t{} \t{} \t{} \t\t {} \t{} \t{} \t{} \t{} \t{} \n'.format(*x))
nnsave.close()

quit()




fig, ax = plt.subplots(figsize=(28,7))

plt.plot(wave_rest, flux_norm,linewidth=1,drawstyle='steps-mid', alpha=0.5, color='grey', label="Observed flux, idx: {}".format(chosen[0]))
plt.plot(wave_rest, sigma_norm,linewidth=1,drawstyle='steps-mid', alpha=0.5, color='grey')


colors=['black','m','blue', 'c']
for j,idx in enumerate(nbrs_idx):
    
    nn_fit = Fits[idx]
    nn_meta_idx = Idxs[idx]
    nn_wave = Wave[idx]

    #plt.plot(wave_rest_array, fit, linewidth=1, drawstyle='steps-mid', alpha=0.4, color=colors[j])
    plt.plot(nn_wave[nint], nn_fit[nint], linewidth=1.5, drawstyle='steps-mid', color=colors[j], label='NN {}, idx: {}'.format(j,nn_meta_idx))
    
#plt.plot(wave_rest_array, chosen_fit, linewidth=1.5, drawstyle='steps-mid', alpha=0.4, color='red')
plt.plot(chosen_wave[nint:b_int], chosen_fit[nint:b_int], linewidth=2, drawstyle='steps-mid', color='red', label='chosen qso, idx: {}'.format(chosen[0]))

plt.xlim(wl_a, wl_b)
plt.ylim(-0.1,max(chosen_fit[nint])+0.4)
plt.legend()
plt.xlabel('$\lambda_{rest}$ [Angstrom]')
plt.ylabel('flux')
plt.title('deltapix2 = {}'.format(dpix))
plt.savefig(dpath + 'plots/{}_{}_{}.png'.format(chosen_info[0], chosen_info[1], chosen_info[2]), format='png')
plt.close()








