import time
import glob
import numpy as np
import imageio as io
from PIL import Image
import matplotlib.pyplot as plt

#togliere i commenti a secnda del necessario

start = time.time()
x = np.loadtxt('datidb.dat', unpack=True)
end = time.time()
print(f'Tempo lettura = {end-start} s')

#infomazioni per il plot
n = int(x[0])
m = int(x[1])
d = int(x[2])
a = int(x[3])

#singolo plot
if a == 0 :
    start = time.time()
    x = x[4:]
    X = np.reshape(x, (n, n))

    for i in range(n):
        for j in range(n):
            if X[i, j]==1: X[i, j]=0
            else: X[i, j]= 1 
        
    plt.figure(1)
    plt.title(f'reticolo {n}x{n}\n numero di punti={m}')
    plt.pcolor(X, cmap='binary')
    
    end = time.time()
    print(f'Tempo per grafico = {end-start} s')
    
    plt.show()
    
    
"""
#tanti plot per l'animazione  
M = int(m/d)
  
if a == 1 :
    start = time.time()
    x = x[4:]
    X = np.reshape(x, (M, n, n))

    for k in range(M):
        for i in range(n):
            for j in range(n):
                if X[k, i, j]==1: X[k, i, j]=0
                else: X[k, i, j]= 1 
        
        fig=plt.figure(k)
        plt.title(f'reticolo {n}x{n}\n numero di punti={m}')
        plt.pcolor(X[k, :, :], cmap='binary')
        plt.savefig(f'/home/francesco/Scrivania/gif/{k}')
        plt.close(fig)
        
    end = time.time()
    print(f'Tempo per salvare {M} grafici = {end-start} s')
"""

"""
#uniorne plot per ottenere la gif      
start = time.time()
    
frames=[]
imgs = sorted(glob.glob('/home/francesco/Scrivania/gif/*.png'))
imgs.sort(key=len)
for i in imgs:
    new_frame = Image.open(i)
    frames.append(new_frame)


frames[0].save('dielectric_breakdown.gif',format='GIF', \
append_images=frames[:],save_all=True,duration=100,loop=0)
    
end = time.time()
print(f"Tempo per produrre l'animazione = {end-start} s")

"""  
