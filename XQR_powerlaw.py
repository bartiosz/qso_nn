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

As = [r2a,r3a,r4a,r5a,r6a,r7a,r8a]
Bs = [r2b,r3b,r4b,r5b,r6b,r7b,r8b]
wl = np.linspace(1200,3000,1000)
mask = np.zeros(len(wl),dtype=bool)
for i,m in enumerate(As):
    mask = np.logical_or(mask, (wl > As[i]) & (wl < Bs[i]))


npath = '/media/bartosz/USB STICK/BOSS_DR14/normed/'
epath = '/media/bartosz/USB STICK/BOSS_DR14/'
