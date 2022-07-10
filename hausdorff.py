import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit as fit

#raggi dei cerchi
rr = [10,20,30,40,50,60]
pp = []

for i in range(1,6):
    ii = str(i)
    file = 'curva'+ii+'.npy'
    curva = np.load(file)
    N = len(curva[0,:])
    c = int((N+1)/2)
    nn = []
    for k, r in enumerate(rr):
        n = 0
        for i in range(c-r,c+r+1):
            for j in range(c-r,c+r+1):
                if (i-c)**2 + (j-c)**2 - r**2 < 0 and curva[i,j] == 0:
                    n = n+1;
        nn.append(n)

    lr = [np.log(x) for x in rr]
    ln = [np.log(x) for x in nn]

    def f(x,a,b):
        return a*x+b

    p, covm = fit(f,lr,ln)
    pp.append(p[0])

d = np.mean(pp)
sd = np.std(pp)

print('%.2f +- %.2f'%(d,sd))

ii = str(5)
file = 'curva'+ii+'.npy'
curva = np.load(file)
path = 'c:/users/gugli/desktop/tesi/simulazioni/plot_curva'+ii+'.pdf'
plt.imsave(fname=path,arr=curva.T,cmap='binary',origin='lower')
