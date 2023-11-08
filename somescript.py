dpath= '/media/bartosz/Volume/BOSS_DR14/data/'

nn = np.loadtxt(dpath + 'nearest_neighbors.txt')

idx = nn[:,0]

nn_idx = nn[:,2:7]


bestN = 110799
bestNNs = [idx[i] for i,nn in enumerate(nn_idx) if bestN in nn]

fig, ax = plt.figure()

specfiles = [fname for fname in glob.glob(dpath + 'fits/dpx25/' + '{}_*.txt'.format(index)) for index in bestNNs]
for fname in specfiles:
	spec = np.loadtxt(fname)
	wl = spec[:,0]
	flux = spec[:,1]

	ax.plot(wl,spec, c='lightgrey', alpha=0.2)


