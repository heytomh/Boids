import sys, argparse
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.spatial.distance import squareform, pdist, cdist
from numpy.linalg import norm


width, height = 640, 480

class Boids:
    """class that represents Boids simulation"""
    def __init__(self, N):
        """initialize the Boid simulation"""
        # initial position and velocities
        self.pos = [width/2.0, height/2.0] + 10*np.random.rand(2*N).reshape(N,2)
        # normalized random velocities
        angles = 2*math.pi*np.random.rand(N)
        self.vel = np.array(list(zip(np.sin(angles), np.cos(angles))))
        self.N = N
        # minimum distance of approach
        self.minDist = 25.0
        # maximum magnitude ofvelocities calculated by "rules"
        self.maxRuleVel = 0.03
        # maximum magnitude of the final velocity
        self.maxVel = 2.0

    def tick(frameNum, pts, beak):
        """update the simulation by one time step."""
        # get pairwise distances
        self.distMatrix = squareform(pdist(self.pos))
        # apply rules:
        self.vel =+ self.applyRules()
        self.limit(self.vel, self.maxVel)
        self.pos += self.vel
        self.applyBC()
        # update data
        pts.set_data(self.pos.reshape(2*self.N)[::2], self.pos.reshape(2 * self.N)[1::2])
        vec = self.pos + 10 * self.vel / self.maxVel
        beak.set_data(vec.reshape(2 * self.N)[::2], vec.reshape(2 * self.N)[1::2])

    def limitVec(self, vec, maxVal):
        """limit the magnitude of the 2D vector"""
        mag = norm(vec)
        if mag > maxVal:
            vec[0], vec[1] = vec[0]*maxVal/mag, vec[1]*maxVal/mag

    def limit(self, X, maxVal):
        """limit the magnitude of 2D vectors in array X to maxValue"""
        for vec in X:
            self.limitVec(vec, maxVal)

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

    def applyRules(self):
        #Rule One - Separation
        D = distMatrix < 25.0
        vel = self.pos*D.sum(axis=1).reshape(self.N, 1) - D.dot(self.pos)
        self.limit(vel, self.maxRuleVel)

        #distance threshold for alignment
        D = distMatrix < 50.0

        #apply rule two - alignment
        vel2 = D.dot(self.vel)
        self.limit(vel2, self.maxRuleVel)
        vel += vel2;

        #apply rule 3 - cohesino
        vel3 = D.dot(self.pos) - self.pos
        self.limit(vel3, self.maxRuleVel)
        vel += vel3

        return vel

    def buttonPress(self, event):
        """event handelr for matplotlib button presses"""
        # left-click to add a boid
        if event.button is 1:
            self.pos = np.concatenate((self.pos, np.array([[event.xdata, event.ydata]])), axis=0)
            # generate random velocity
            angles = 2 * math.pi * np.random.rand(1)
            v = np.array(list(zip(np.sin(angles), np.cos(angles))))
            self.vel = np.concatenate((self.vel, v), axis=0)
            self.N += 1
        # right-click to scatter boids
        elif event.button is 3:
            # add scattering velocity
            self.vel += 0.1 * (self.pos - np.array([[event.xdata, event.ydata]]))

def tick(frameNum, pts, beak, boids):
    boids.tick(frameNum, pts, boids)
    return pts, beak

# main() function
def main():
    # use sys.argv if needed
    print('starting boids...')

    parser = argparse.ArgumentParser(description="Implementing Craig Reynolds's Boids...")

    # add arguments
    parser.add_argument('--num-boids', dest='N', required=False)
    args = parser.parse_args()

    # set the initial number of boids
    N = 100
    if args.N:
        n = int(args.N)

    # create boids
    boids = Boids(N)

    # set up plot
    fig = plt.figure()
    ax = plt.axes(xlim=(0,width), ylim=(0,height))

    pts, = ax.plot([],[],markersize=10, c='k', marker='o', ls='None')
    beak, = ax.plot([],[],markersize=4, c='r', marker='o', ls='None')
    anim = animation.FuncAnimation(fig, tick, fargs=(pts, beak, boids), interval=50)

    # add a "button press" event handler
    cid = fig.canvas.mpl_connect('button_press_event', buttonPress)

    plt.show()

# call main
if __name__ == '__main__':
    main()