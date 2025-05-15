import numpy as np
import matplotlib.pyplot as plt

import sys, getopt

from astropy import units as u
from astropy.coordinates import SkyCoord

from specdb import query_catalog as spqcat
from specdb import utils as spdbu
from specdb.specdb import SpecDB, IgmSpec
from specdb import specdb as sdbsdb



def degtodegms(c):
    deg,cc=divmod(c,1)
    m,s=divmod(cc*60,1)
    return '{:.0f}° {:.0f}\' {:.1f}\'\''.format(deg,m,s)

def degtohms(c):
    h,cc=divmod(c/360*24,1)
    m,s=divmod(cc*60,1)
    return '{:.0f}h {:.0f}\' {:.1f}\'\''.format(h,m,s)


def main(argv):
    ra, dec, z = '-','-','-'
    opts, args = getopt.getopt(argv,"hi:s:",["index=","save="])
    for opt, arg in opts:
        if opt == '-h':
            print ('plot_spectrum.py -i <object index> -s <y/n>')
            sys.exit()
        elif opt in ("-i", "--index"):
            idx = int(arg)
            igmsp = SpecDB('/media/bartosz/Volume/igmspec_data/DB/IGMspec_DB_v03.1.hdf5')
            meta = igmsp['BOSS_DR14'].meta
            obj = meta[idx]
            ra,dec,z = obj[0],obj[1],obj[6]
            spec_file=obj[-1]


            coords = SkyCoord(ra,dec,unit='deg')

            meta_sing = igmsp.meta_from_position(coords, 1*u.arcsec)
            spec = igmsp.spectra_from_meta(meta_sing)
            
            fig, ax = plt.subplots(figsize=(20,7))
            
            ax.plot(spec.wavelength,spec.sig, lw=1,drawstyle='steps-mid', color='grey',label='errors')
            ax.plot(spec.wavelength,spec.flux, lw=2,c='black',drawstyle='steps-mid',label='observed flux')

            info = '\n'.join((
                r'RA = {}'.format(degtohms(ra)),
                r'DEC = {}'.format(degtodegms(dec)),
                r'$z$ = {0:.3f}'.format(z)))

            props = dict(boxstyle='square', facecolor=None, alpha=0.5)
            ax.text(0.75, 0.97, info, transform=ax.transAxes, fontsize=12,
                verticalalignment='top', bbox=props)

            plt.legend(fontsize=12)
            plt.xlabel('$\lambda (\AA)$')
            plt.ylabel('flux')

        elif opt in ('-s', '--save'):
            if arg == 'y':
                plt.savefig('plots/{}.png'.format(spec_file[:-5]),format='png')

        
    print ('RA=', ra, ' ; DEC=', dec, ' ; Z=', z, 'spec_file=', spec_file)
    plt.show()
    plt.close()

if __name__ == "__main__":
    main(sys.argv[1:])
