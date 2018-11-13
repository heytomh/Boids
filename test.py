# vec = self.pos + 10*self.vel/self.maxVel
# beak.set_data(vec.reshape(2*self.N)[::2], ved.reshape(2*self.N)[1::2])

def test2(pos, radius):
    distMatrix = squareform(pdist(pos))
    D = distMatrix < radius
    vel = pos*D.sum(axis=1).reshape(N,1) - D.dot(pos)
    return vel