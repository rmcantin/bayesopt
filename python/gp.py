import numpy as np

class GaussianProcess:
    def __init__():
        self.x = np.array([])
        self.y = np.array([])
    
    def correlationFunction(x1,x2,theta):
        xdiff = x1 - x2

        for xd in xdiff:
            dif = abs(xd) / theta
            sum += dif
            prod *= 1+dif

        return prod*exp(-sum)

    def computeGPparams():
        uno = ones(N,1)

        uInvR = dot(uno.T,invR)
        eta = dot(uInvR,uno) + 1/delta

        mu = dot(uInvR,normY) / eta;
        yRy = dot(dot(normY.T,invR),normY)
        sig2 = (b + yRy - mu*mu/eta) / (alpha+N+2)

        return mu,sig2

    def prediction(query):
        colR = zeros(nPoints,1)
        
        for rows in self.x:
            
    
