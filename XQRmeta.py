import numpy as np

mpath = '/media/bartosz/Volume/XQR30/data/'

meta = np.loadtxt(mpath+'meta.txt', dtype='str')

names = meta[:,0]
Z = [float(z) for z in meta[:,1]]

wl1 = 13450
wl2 = 14250

wl1r = [wl1/(1+z) for z in Z]
wl2r = [wl2/(1+z) for z in Z]

meta_file = mpath +'meta.txt'
with open(meta_file,'w') as msave:
    for x in zip(names,Z,wl1r,wl2r):
        msave.write('{} \t {} \t {} \t {} \n'.format(*x))
msave.close()
