import numpy as np

mpath = '/media/bartosz/Volume/highz_data/'

meta = np.loadtxt(mpath+'meta_data.txt', dtype='str')

names = meta[:,0]
Z = [float(z) for z in meta[:,1]]
SNR = [float(snr) for snr in meta[:,2]]

wl1 = 13450
wl2 = 14250

wl1r = [wl1/(1+z) for z in Z]
wl2r = [wl2/(1+z) for z in Z]

meta_file = mpath +'meta_data_v2.txt'
with open(meta_file,'w') as msave:
    for x in zip(names,Z,SNR,wl1r,wl2r):
        msave.write('{} \t {} \t {} \t {} \t {} \n'.format(*x))
msave.close()
