import numpy as np
import glob
import matplotlib.pyplot as plt
from specdb.specdb import SpecDB, IgmSpec

dpath = '/media/bartosz/Volume/BOSS_DR14/data/'
fpath = dpath + 'fits/dpx25/'
npath = dpath + 'normed/dpx25/'

def plot_data(idx,fig,ax,c,lw):

    filename = glob.glob(npath + '{}_*_0_norm.txt'.format(idx))[0]
    normed = np.loadtxt(filename)

    norm_wl = normed[:,0]
    norm_flux = normed[:,1]
    norm_sig = normed[:,2]

    ax.plot(norm_wl,norm_sig,linewidth=0.8*lw,drawstyle='steps-mid',color=c,alpha=0.5,label='errors')
    ax.plot(norm_wl,norm_flux,linewidth=lw, drawstyle='steps-mid', color=c,alpha=0.8,label='normed observation')
    return 
    

def plot_fit(idx,fig,ax,c,lw):

    filename = glob.glob(fpath + '{}_*_0_dpx25.txt'.format(idx))[0]
    fit = np.loadtxt(filename)

    fit_wl = fit[:,0]
    fit_flux = fit[:,1]

    mask1 = fit_wl > 1250
    mask2 = fit_wl < 2250
    mask = np.logical_and(mask1,mask2)


    #for filename in glob.glob(npath + '{}_*_0_norm.txt'.format(idx)):
    #    normed = np.loadtxt(filename)

    #norm_wl = normed[:,0]
    #norm_flux = normed[:,1]
    #norm_sig = normed[:,2]


    #fig = plt.figure(figsize=(20,7))
    #ax = fig.add_subplot(1,1,1)

    #ax.plot(norm_wl,norm_sig,linewidth=0.8,drawstyle='steps-mid',color='grey',label='errors')
    #ax.plot(norm_wl,norm_flux,linewidth=1, drawstyle='steps-mid', color='black',label='normed observation')
    ax.plot(fit_wl,fit_flux,linewidth=lw, drawstyle='steps-mid', color=c,label='spline fit')
    #ax.axvspan(min(min(norm_wl),min(fit_wl)),1250,alpha=0.2,facecolor='grey')
    #ax.axvspan(2250,max(max(norm_wl),max(fit_wl)),alpha=0.2,facecolor='grey')
    #ax.set_xlim(left=min(min(norm_wl),min(fit_wl)),right=max(max(norm_wl),max(fit_wl)))
    return max(fit_flux[mask])


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
