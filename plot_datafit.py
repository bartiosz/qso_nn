import matplotlib.pyplot as plt
import numpy as np
from specdb import query_catalog as spqcat
#from specdb import interface_db as spgidb
from specdb import utils as spdbu
from specdb.specdb import SpecDB, IgmSpec
from specdb import specdb as sdbsdb
import importlib
import pickle
import glob
import sys, getopt


def main(argv):
    dpath = '/media/bartosz/Volume/BOSS_DR14/data/'
    fpath = dpath + 'fits/dpx25/'
    npath = dpath + 'normed/dpx25/'
    ra, dec, z, spec_file = '-','-','-','-'
    opts, args = getopt.getopt(argv,"hi:s:",["index=","save="])
    for opt, arg in opts:
        if opt == '-h':
            print ('plot_spectrum.py -i <object index> -s <y/n>')
            sys.exit()
        elif opt in ("-i", "--index"):
            idx = int(arg)
            #igmsp = SpecDB('/media/bartosz/Volume/igmspec_data/DB/IGMspec_DB_v03.1.hdf5')
            #meta = igmsp['BOSS_DR14'].meta
            #obj = meta[idx]
            #ra,dec,z = obj[0],obj[1],obj[6]
            #spec_file=obj[-1]


            for filename in glob.glob(fpath + '{}_*_0_dpx25.txt'.format(idx)):
                fit = np.loadtxt(filename)

            fit_wl = fit[:,0]
            fit_flux = fit[:,1]


            for filename in glob.glob(npath + '{}_*_0_norm.txt'.format(idx)):
                normed = np.loadtxt(filename)

            norm_wl = normed[:,0]
            norm_flux = normed[:,1]
            norm_sig = normed[:,2]

      
            fig, ax = plt.subplots(figsize=(20,7))
            fig = plt.figure(figsize=(20,7))
            ax = fig.add_subplot(1,1,1)

            ax.plot(norm_wl,norm_sig,linewidth=1,drawstyle='steps-mid',color='grey',label='errors')
            ax.plot(norm_wl,norm_flux,linewidth=2, drawstyle='steps-mid', color='black',label='normed observation')
            ax.plot(fit_wl,fit_flux,linewidth=2, drawstyle='steps-mid', color='red',label='spline fit')
            #plt.show()

        elif opt in ('-s', '--save'):
            if arg == 'y':
                plt.savefig('plots/{}_fit.png'.format(specfile[:-5]),format='png')
            else:
                fname = arg
                plt.savefig('plots/{}_{}.png'.format(idx,fname),format='png')

    return fig
        
    #print ('RA=', ra, ' ; DEC=', dec, ' ; Z=', z, 'spec_file=', spec_file)
    #plt.show()
    plt.close()

if __name__ == "__main__":
    main(sys.argv[1:])

