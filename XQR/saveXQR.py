import glob
import matplotlib.pyplot as plt
from astropy.io import fits
from pathlib import Path
import numpy as np
import itertools


xpath = '/media/bartosz/Volume/XQR30/'
xqpath = xpath + 'XQR30_latest/'
spath = xpath + 'data/spec/'


flist = [f for f in glob.glob(xqpath + '*/*_rebinned_50kms_spec.fits')]

pred = np.loadtxt(xpath + 'pred_masterfile_final_bugfix',dtype='str')
Z = pred[:,1]
qnames = [p.split('.')[0].split('_')[0] for p in pred[:,2]]

msave = open(xpath + 'data/meta.txt','w')

for f in flist:
    ffile_name = Path(f).stem
    ffile_info = ffile_name.split('_')
    qname = ffile_info[0]
    print(qname)
    if qname in qnames:
        with fits.open(f) as hdul:
            data = hdul[1].data

        wl = np.array([d[0] for d in data])*10    #convert from nm to Angstrom
        flux = np.array([d[3] for d in data])
        sig = np.array([d[4] for d in data])

        mask1 = sig < np.median(sig)*10
        mask2 = sig > -np.median(sig)*10
        mask = np.logical_and(mask1,mask2)
       
        #wl = wl[mask]
        #flux = flux[mask]
        #sig = sig[mask]

        spec_save = spath + '{}.txt'.format(ffile_name)
        with open(spec_save,'w') as ssave:
            for x in zip(wl,flux,sig):
                ssave.write('{} \t {} \t {} \n'.format(*x))
        ssave.close()

        qidx = qnames.index(qname)
        z = Z[qidx]
        msave.write('{} \t {} \n'.format(qname,z))        

msave.close()
