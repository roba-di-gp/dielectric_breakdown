import numpy as np
import matplotlib.pyplot as plt
from numba import jit
from time import time
start = time()

@jit
def solve(n, iters, ground, seed):
    #inizializzo il potenziale con le condizioni al contorno
    phi = np.ones([n,n])
    for i in range(1,n-1):
        for j in range(1,n-1):
            phi[i,j] = seed[i,j]
            if ground[i,j] == 0:
                phi[i,j] = 0

    #risolvo l'equazione di laplace
    for k in range(iters):
        for i in range(1,n-1):
            for j in range(1,n-1):
                if ground[i,j] == 0:
                    phi[i,j] = 0;
                else:
                    phi[i,j] = (phi[i,j+1]+phi[i,j-1]+phi[i+1,j]+phi[i-1,j])/4
    return phi
