import numpy as np
from math import factorial
from numpy.matlib import dot,ones,array
from scipy.stats import norm
#norm.pdf(x, loc=0, scale=1) and norm.cdf(x, loc=0, scale=1)

class Kernel:
    def __init__(self):
        self.theta = 0.21
        self.p = 1.6

    def mattern(self,x1,x2):
        theta = self.theta
        xdiff = x1 - x2

        acsum = 0
        acprod = 1

        for xd in xdiff:
            dif = abs(xd) / theta
            acsum += dif
            acprod *= 1+dif

        return acprod*exp(-acsum)

    
class GaussianProcess:
    def __init__(self):
        self.x = array([])
        self.y = array([])
        self.ny = array([])
        self.min = 0
        self.max = 0
        
        self.prior_alpha = 1.0
        self.prior_beta = 0.1
        self.prior_delta = 10.0

        self.gp_noise = 1e-4

    
    def correlationFunction(self,x1,x2):
        k = kernel()
        return k.mattern(x1,x2)

    def addNewPoint(self,newX,newY):
        self.x = np.r_[self.x,newX]
        self.y = np.r_[self.y,newY]

        minIndex = self.y.argmin()
        maxIndex = self.y.argmax()

        normalize = False

        if minIndex != self.min:
            self.min = minIndex
            normalize = True

        if maxIndex != self.max:
            self.max = maxIndex
            normalize = True

        if normalize:
            normalizeResponseData()

    def normalizeResponseData(self):
        valMin = self.y[self.min]
        valMax = self.y[self.max]
        
        self.ny = (self.y - valMin) / (valMax-valMin)

    def mahalanobisDistance(self,m,invC):
        return dot(dot(m,invC),m)

    def computeGPparams(self):
        N = self.y.shape
        uno = ones(N,1)  #TODO: Generalize for other mean functions
        alpha = self.prior_alpha
        beta  = self.prior_beta
        delta = self.prior_delta

        uInvR = dot(uno,invR)
        eta = dot(uInvR,uno) + 1/delta
        yRy = mahalanobisDistance(self.ny,invR)
        
        mu = dot(uInvR,self.ny) / eta;
        sig2 = ( beta + yRy - mu*mu/eta ) / (alpha+N+2)

        return mu,sig2


    def prediction(self,query,mu,sig2):
        uno = ones(N,1)  #TODO: Generalize for other mean functions
        r = array([])
        for x in dataX:
            r = np.r_[r,correlationFunction(x,query,theta)]

        rInvR = dot(r,invR)
        uInvR = dot(uno,invR)

        uInvRr = dot(uInvR,r)
        rInvRr = dot(rInvR,r)

        ypred  = mu + rInvR * (normY - uno * mu)
        spred2 = sig2 + (1 - rInvRr + (1 - uInvRr)**2 / eta )

        return ypred,spred2

    def run(self):
        print self.gp_noise


class SKO:
    def __init__(self):
        self.params_SKO_method = 'EI'
        self.params_SKO_maxiter = 300
        self.params_SKO_stop = 1e-4
        self.params_SKO_maxdims = 20

        self.params_inner_maxeval = 1000
        self.params_inner_maxiter = 300

        self.params_lhs_evalsdim = 10
        self.params_lhs_maxevals = 100
        


if __name__ == '__main__':
    myopt = GaussianProcess()
    myopt.run()
    
