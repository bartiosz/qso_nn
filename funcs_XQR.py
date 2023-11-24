import numpy as np
import glob
import matplotlib.pyplot as plt
from specdb.specdb import SpecDB, IgmSpec

dpath = '/media/bartosz/Volume/XQR30/data/'
fpath = dpath + 'fits/'
npath = dpath + 'normed/'

def plot_xdata(idx,fig,ax,c,lw):
    filename = glob.glob(npath + '{}_*_norm.txt'.format(idx))[0]
    normed = np.loadtxt(filename)

    norm_wl = normed[:,0]
    norm_flux = normed[:,1]
    norm_sig = normed[:,2]
    ax.plot(norm_wl,norm_sig,linewidth=lw*0.8,drawstyle='steps-mid',color=c,alpha=0.5,label='errors')
    ax.plot(norm_wl,norm_flux,linewidth=lw, drawstyle='steps-mid', color=c,alpha=0.8,label='normed observation')
    
    return

def plot_xfit(idx,fig,ax,c,lw):

    filename = glob.glob(fpath + '{}_*_dpx34.txt'.format(idx))[0]
    fit = np.loadtxt(filename)

    fit_wl = fit[:,0]
    fit_flux = fit[:,1]

    mask1 = fit_wl > 1250
    mask2 = fit_wl < 2250
    mask = np.logical_and(mask1,mask2)

    ax.plot(fit_wl,fit_flux,linewidth=lw, drawstyle='steps-mid', color=c,label='spline fit')
    #ax.axvspan(min(min(norm_wl),min(fit_wl)),1250,alpha=0.2,facecolor='grey')
    #ax.axvspan(2250,max(max(norm_wl),max(fit_wl)),alpha=0.2,facecolor='grey')
    #ax.set_xlim(left=min(min(norm_wl),min(fit_wl)),right=max(max(norm_wl),max(fit_wl)))
    return max(fit_flux[mask])

