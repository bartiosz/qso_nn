import numpy as np
import glob
import matplotlib.pyplot as plt
from specdb.specdb import SpecDB, IgmSpec

dpath = '/media/bartosz/Volume/BOSS_DR14/data/'
fpath = dpath + 'fits/dpx25/'
npath = dpath + 'normed/dpx25/'

def plot_fit(idx):

    for filename in glob.glob(fpath + '{}_*_0_dpx25.txt'.format(idx)):
        fit = np.loadtxt(filename)

    fit_wl = fit[:,0]
    fit_flux = fit[:,1]


    for filename in glob.glob(npath + '{}_*_0_norm.txt'.format(idx)):
        normed = np.loadtxt(filename)

    norm_wl = normed[:,0]
    norm_flux = normed[:,1]
    norm_sig = normed[:,2]


    fig = plt.figure(figsize=(20,7))
    ax = fig.add_subplot(1,1,1)

    ax.plot(norm_wl,norm_sig,linewidth=1,drawstyle='steps-mid',color='grey',label='errors')
    ax.plot(norm_wl,norm_flux,linewidth=2, drawstyle='steps-mid', color='black',label='normed observation')
    ax.plot(fit_wl,fit_flux,linewidth=2, drawstyle='steps-mid', color='red',label='spline fit')
    ax.axvspan(min(min(norm_wl),min(fit_wl)),1250,alpha=0.2,facecolor='grey')
    ax.axvspan(2250,max(max(norm_wl),max(fit_wl)),alpha=0.2,facecolor='grey')
    ax.set_xlim(left=min(min(norm_wl),min(fit_wl)),right=max(max(norm_wl),max(fit_wl)))
    return fig, ax


def get_meta():
    igmsp = SpecDB('/media/bartosz/Volume/igmspec_data/DB/IGMspec_DB_v03.1.hdf5')
    meta = igmsp['BOSS_DR14'].meta
    return meta

def get_z(meta,idx):
    obj = meta[idx]
    z = obj[6]
    return z

def get_BI(meta,idx):
    obj = meta[idx]
    BI = obj[27]
    return BI
