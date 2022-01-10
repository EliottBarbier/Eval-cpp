import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation


tab=np.genfromtxt("build/Euler_explicit.txt") #Ou implicite
print(tab.shape)
temps=tab[:,0]
T=[]

for k in range(tab.shape[0]):
    T.append(tab[k,1:])

Nx = 50 #On reporte le nombre de points que l'on a pris : pas obligatoire, pourrait se trouver dans les températures.
x_max=1 #Ainsi que le x_max 

deltax=x_max/Nx

X=[k*deltax for k in range(Nx)]

fig = plt.figure() # initialise la figure
line, = plt.plot(X, T[0]) 
plt.xlim(0,1)


def animate(i):
    line.set_data(X, T[i])
    #line.set_label(f"t={temps[i]}")
    #plt.legend()
    return line,
    
ani = animation.FuncAnimation(fig, animate, frames=len(T), blit=True, interval=100 repeat=True)
plt.show()