import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

tab=np.genfromtxt("build/Euler_explicite.txt")


print(tab.shape) #Vérification
temps=tab[:,0] #Le tableau des temps
T=[] #La liste où l'on va mettre les vecteurs colonnes des températures à tout instant.

for k in range(tab.shape[0]):
    T.append(tab[k,1:])

Nx = tab.shape[1]-1
x_max=1. #Ne pas oublier de la changer si on le change dans le main.cpp
deltax=x_max/Nx

X=[k*deltax for k in range(Nx)] #La liste des x_i


fig = plt.figure() # On initialise la figure
ax = fig.add_subplot(111, autoscale_on=True,
                     xlim=(0, 1.), ylim=(-2, 2))
line, = ax.plot([], [])
time_text = ax.text(0.02, 0.95, '', transform=ax.transAxes)                    

line.set_data(X, T[0])
time_text.set_text('time = 0')


def animate(i): #La fonction qui permet de modifier le plot et le temps qui est affiché, i correspond au numéro de l'image.
    line.set_data(X, T[i])
    time_text.set_text(f'time = {temps[i]:.4f}')
    return line, time_text
ax.set_xlabel("Position x")
ax.set_ylabel("Température")
ax.legend()
anim = animation.FuncAnimation(fig, animate, frames=len(T), blit=True, interval=100, repeat=True) #Pour animer le graphe de la température

plt.show()