import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation


width, height = 640, 480

pos = [width/2.0, height/2.0] + 10*np.random.rand(2*N).reshape(N,2)
angles = 2*math.pi*np.random.rand(N)
vel = np.array(list(zip(np.sin(angles), np.cos(angles))))

def applyBC(self):
    """apply boundry conditions"""
    deltaR = 2.0
    for coord in self.pos:
        if coord[0] > width + deltaR:
            coord[0] = - deltaR
        if coord[0] < - deltaR:
            coord[0] = width + deltaR
        if coord[1] > height + deltaR:
            coord[1] = - deltaR
        if coord[1] < - deltaR:
            coord[1] = height + deltaR

fig = plt.figure()
ax = plt.axes(xlim=(0,width), ylim=(0,height))

pts, = ax.plot([],[],markersize=10, c='k', marker='o', ls='None')
beak, = ax.plot([],[],markersize=4, c='r', marker='o', ls='None')
anim = animation.FuncAnimation(fig, tick, fargs=(pts, beak, boids), interval=50)

