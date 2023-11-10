import matplotlib.pyplot as plt
from matplotlib import colors
import numpy as np
import importlib
import glob
import sys, getopt
import os
import glob
from matplotlib.backends.backend_pdf import PdfPages
from funcs import plot_fit, get_meta, get_z, get_BI

def take_second(elem):
    return elem[1]


dpath= '/media/bartosz/Volume/BOSS_DR14/data/'

nn = np.loadtxt(dpath + 'nearest_neighbors.txt')

idx = nn[:,0]

nn_idx = nn[:,2:7]
nn_d = nn[:,8:13]

d_mean = np.array([np.mean(d) for d in nn_d])

idx_dmean = sorted([[ii,d_mean[i]] for i,ii in enumerate(idx)],key=take_second)

top500 = idx_dmean[:]


meta = get_meta()
pp = PdfPages('plots/qso0_99.pdf')
for i,qso in enumerate(top500):
    qidx = int(qso[0])
    qdmean = qso[1]
    z, BI = get_z(meta,qidx), get_BI(meta,qidx)

    fig, ax = plot_fit(qidx)
    ax.set_xlabel('$\lambda_{rest} (\AA)$')
    ax.set_ylabel('flux')
    ax.legend()

    info = '\n'.join((
        r'id = {}'.format(qidx),
        r'$d_{{mean}}$ = {0:.3f}'.format(qdmean),
        r'BI = {0:.1f}'.format(BI),
        r'$z$ = {0:.3f}'.format(z)))


    
    fc = colors.to_rgba('white')
    ec = colors.to_rgba('grey')        
    fc = fc[:-1] + (0.75,)
    ec = ec[:-1] + (0.15,)
    props = dict(boxstyle='square', facecolor=fc, edgecolor=ec)
    ax.text(0.75, 0.975, info, transform=ax.transAxes, fontsize=12,
        verticalalignment='top', bbox=props)

    #plt.show()
    #quit()
    pp.savefig(fig)

    if (i+1) % 100 == 0:
        pp.close()
        pp = PdfPages('plots/qso{}_{}.pdf'.format(i+1,i+100))
    
    print(i)

pp.close()

