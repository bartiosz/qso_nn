import sys, getopt
from specdb.specdb import SpecDB, IgmSpec

def main(argv):
    ra, dec, z = '-','-','-'
    opts, args = getopt.getopt(argv,"hi:",["index="])
    for opt, arg in opts:
        if opt == '-h':
            print ('object_info.py -i <object index> ')
            sys.exit()
        elif opt in ("-i", "--index"):
            idx = int(arg)
            igmsp = SpecDB('/media/bartosz/Volume/igmspec_data/DB/IGMspec_DB_v03.1.hdf5')
            meta = igmsp['BOSS_DR14'].meta
            obj = meta[idx]
            ra,dec,z = obj[0],obj[1],obj[6]
            spec_file=obj[-1]

        
    print ('RA=', ra, ' ; DEC=', dec, ' ; Z=', z, 'spec_file=', spec_file)

if __name__ == "__main__":
    main(sys.argv[1:])
