import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from pathlib import Path
import glob
import itertools

from scipy.optimize import curve_fit

def power_law(x,alpha,beta):
    return alpha*x**beta

r1a, r1b = 1150,1170
r2a, r2b = 1275,1290
r3a, r3b = 1350,1360
r4a, r4b = 1445,1465
r5a, r5b = 1690,1705
r6a, r6b = 1770,1810
r7a, r7b = 1970,2400
r8a, r8b = 2480,2675
r9a, r9b = 2925,3400

As = [r2a,r3a,r4a,r5a,r6a,r7a,r8a,r9a]
Bs = [r2b,r3b,r4b,r5b,r6b,r7b,r8b,r9b]

dpath = '/media/bartosz/USB STICK/highz_data/'
npath = dpath + 'fits/'

# mask water absorption
wl1 = 13450
wl2 = 14250


spec_id = []
alpha,beta,alpha_sig,beta_sig = [],[],[],[]
rsq = []

meta = np.loadtxt(dpath + 'meta_data_v2.txt', dtype='str')
meta_names = meta[:,0]
meta_z = [float(z) for z in meta[:,1]]


pp = PdfPages('plots/XQR_powerlaw_spline_test3.pdf')

for i,f in enumerate(glob.glob(npath + '*.txt')):    

    norm = np.loadtxt(f)
#    wl,flux,sig = norm[:,0],norm[:,1],norm[:,2]
    wl,flux = norm[:,0],norm[:,1]
    
    file_name = Path(f).stem
    file_info = file_name.split('_')
    
    spec_id.append(file_info[0])

    meta_id = np.where(meta_names==file_info[0])[0][0]
    qso_z = meta_z[meta_id]
    
    mask = np.zeros(len(wl),dtype=bool)

    if 0 in flux[wl>r9a]:
        N = len(As)-1
    else: N = len(As)

    for j in range(N):
        mask = np.logical_or(mask, (wl > As[j]) & (wl < Bs[j]))
    water_mask = np.logical_or(wl<wl1/(1+qso_z),wl>wl2/(1+qso_z))
    mask = np.logical_and(mask,water_mask)
#    mask = np.logical_and(mask, sig>0)

#    # mask uncertainties >2sigma
#    sig_mean = np.mean(sig[mask])
#    sig_stdev = np.std(sig[mask])
#    mask = np.logical_and(mask,sig<sig_mean+3*sig_stdev)

    try:    
#        popt,pcov = curve_fit(power_law,wl[mask],flux[mask],sigma=sig[mask],absolute_sigma=True,maxfev=1000)
        popt,pcov = curve_fit(power_law,wl[mask],flux[mask],maxfev=2000)
        
    except RuntimeError:
        print(i, 'RuntimeError: 1st iteration')
        continue
        
    diff = abs(flux - power_law(wl,*popt))
    diff_mean = np.mean(diff[mask])
    diff_stdev = np.std(diff[mask])
    mask2 = np.logical_and(mask,diff<diff_mean+3*diff_stdev)

    try:
        #popt,pcov = curve_fit(power_law,wl[mask2],flux[mask2],sigma=sig[mask2],absolute_sigma=True,maxfev=1000)
        popt,pcov = curve_fit(power_law,wl[mask2],flux[mask2],bounds=([popt[0],-np.inf],[np.inf,popt[1]]),maxfev=2000)

    except RuntimeError:
        print(i, 'RuntimeError: 2nd iteration')
        continue

    res = np.sum( (flux[mask]-power_law(wl[mask],*popt))**2 )
    tot = np.sum( (flux[mask]-np.mean(flux[mask]))**2 )
    r_squared = 1 - res/tot

    
    alpha.append(popt[0])
    beta.append(popt[1])
    alpha_sig.append(np.sqrt(pcov[0][0]))
    beta_sig.append(np.sqrt(pcov[1][1]))

    rsq.append(r_squared)


#    median = np.median(flux)
#    stdev = np.std(flux)
#    pmask = abs(flux)<median+5*stdev
    
    fig = plt.figure()
    fig.set_figheight(4.8)
    fig.set_figwidth(12)
    plt.title(file_info[0])
#    plt.plot(wl[pmask],flux[pmask],drawstyle='steps-mid',lw=0.5)
#    plt.plot(wl[pmask],power_law(wl[pmask],*popt),alpha=0.7,lw=0.5)
#    plt.plot(wl[pmask],flux[pmask]-power_law(wl[pmask],*popt),drawstyle='steps-mid',alpha=0.7,lw=0.5)
#    plt.plot(wl[pmask],sig[pmask]-1,drawstyle='steps-mid',alpha=0.5,lw=0.5)
    plt.plot(wl,flux,drawstyle='steps-mid',lw=0.5)
    plt.plot(wl,power_law(wl,*popt),alpha=0.7,lw=0.5)
    plt.plot(wl,flux-power_law(wl,*popt),drawstyle='steps-mid',alpha=0.7,lw=0.5)
    for j,a in enumerate(As):
        plt.axvspan(As[j],Bs[j],color='lightgrey')
    plt.xlabel('wavelength')
    pp.savefig(fig)
    plt.close()
    
    print(i, qso_z)

pp.close()

pl_save='power_law_fits_spline_test3'
with open(dpath + pl_save + '.txt','w') as plsave:
    for x in itertools.zip_longest(spec_id,alpha,alpha_sig,beta,beta_sig,rsq):
        plsave.write('{} \t {} \t {} \t {} \t {} \t {} \n'.format(*x))
plsave.close()
