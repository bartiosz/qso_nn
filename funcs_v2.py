import numpy as np

dpath = '/media/bartosz/Volume/BOSS_DR14/'
meta_file = dpath + 'meta_data.txt'
dist_file = dpath + 'nearest_neighbors.txt'

def get_z(qso_idx):
    meta = np.loadtxt(meta_file, dtype='str')
    qso_IDs = [int(ID) for ID in meta[:,0]]
    list_index = qso_IDs.index(qso_idx)
    zs = meta[:,1]
    z = zs[list_index]
    return float(z)


def get_SNR(qso_idx):
    meta = np.loadtxt(meta_file, dtype='str')
    qso_IDs = [int(ID) for ID in meta[:,0]]
    list_index = qso_IDs.index(qso_idx)
    SNRs = meta[:,2]
    SNR = SNRs[list_index]
    return float(SNR)

def get_meta(qso_idx):
    meta = np.loadtxt(meta_file, dtype='str')
    qso_IDs = [int(ID) for ID in meta[:,0]]
    list_index = qso_IDs.index(qso_idx)
    zs = meta[:,1]
    z = zs[list_index]
    SNRs = meta[:,2]
    SNR = SNRs[list_index]
    return float(z), float(SNR)


def get_dmean(qso_idx):
    meta = np.loadtxt(dist_file)
    qso_IDs = [int(ID) for ID in meta[:,0]]
    list_index = qso_IDs.index(qso_idx)
    dist = (meta[list_index])[8:]
    dmean = np.mean(dist)
    return float(dmean)
