import numpy as np
import matplotlib.pyplot as plt
from numba import jit
from time import time
from numpy.random import choice
from laplace import solve
start = time()

eta = 1 #parametro di controllo del campo
n = 30 #dimensione della matrice
points = 150 #punti della curva di rottura oltre al seme
c = int((n+1)/2) #coordinata del centro

#seme della curva di rottura
ground = np.ones([n,n])
ground[c,c] = 0

#calcolo del potenziale per la prima volta
seed = 0.5*np.ones([n,n])
phi = solve(n, 3000, ground, seed)

for l in range(points):
    #liste dei punti prossimi alla curva, delle loro coordinate e del potenziale
    prox = []; xx = []; yy = []; pot = []
    #elenco i punti prossimi e memorizzo le loro coordinate e il potenziale
    k = 0
    for i in range(1,n-1):
        for j in range(1,n-1):
            if phi[i,j] != 0:
                if ground[i-1,j] == 0 or  ground[i+1,j] == 0 or ground[i,j-1] == 0 or ground[i,j+1] == 0:
                    #conto tutti i legami possibili
                    c = 0
                    if ground[i+1,j] == 0:
                        c = c+1;
                    if ground[i-1,j] == 0:
                        c = c+1;
                    if ground[i,j+1] == 0:
                        c = c+1;
                    if ground[i,j-1] == 0:
                        c = c+1;

                    prox.append(k)
                    xx.append(i)
                    yy.append(j)
                    pot.append(c**(-eta)*phi[i,j])
                    k = k+1;

    if l % 100 == 0:
        print('%.2f'%(l/points))
    #scelgo un punto con la distribuzione pesata
    pp = [x**eta for x in pot]
    norm = sum(pp)
    pp = pp/norm
    ch = choice(prox,p = pp)

    #aggiungo il punto alla curva di rottura
    ground[xx[ch],yy[ch]] = 0

    #ricalcolo il potenziale a partire da quello precedente
    phi = solve(n, 100, ground, phi)

print('\n')
print(time() - start)
print('\n')

#plot
fig0 = plt.figure(0)
plt.imshow(ground.T,cmap='binary',origin='lower')

plt.imsave(fname=r'c:/users/gugli/desktop/tesi/simulazioni/esempio.pdf',arr=ground.T,cmap='binary',origin='lower')
plt.imsave(fname=r'c:/users/gugli/desktop/tesi/simulazioni/esempio_phi.pdf',arr=phi.T,cmap='hot',origin='lower')

np.save('c:/users/gugli/desktop/tesi/simulazioni/esempio.npy',ground)

plt.show()

