import numpy as np
import pylab as pl
import scipy as sp
from math import factorial
from numpy.matlib import dot,ones,array,eye
from scipy.stats import norm
#norm.pdf(x, loc=0, scale=1) and norm.cdf(x, loc=0, scale=1)

class Kernel:
    def __init__(self,params=None):
        if params == None:
            self.set_default_params();

    def set_default_params(self):
        self.params = {'theta':0.10,
                       'p':1.6}

    def mattern(self,x1,x2):
        theta = self.params['theta']   # sqrt(3) / l
        r = abs(x1 - x2)
        dif = r / theta
        acsum = dif.sum()
        acprod = np.prod(1+dif)

        return acprod*np.exp(-acsum)

    
class GaussianProcess:
    def __init__(self):
        self.x = array([])
        self.y = array([])
        self.ny = array([])
        self.min = 0
        self.max = 0
        
        self.prior_alpha = 3.0
        self.prior_beta = 0.1
        self.prior_delta = 10000.0

        self.obs_noise = 1e-8
        self.k = Kernel()
    
    def correlationFunction(self,x1,x2):
        return self.k.mattern(x1,x2)

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
            self.normalizeResponseData()

    def normalizeResponseData(self):
#        valMin = self.y[self.min]
#        valMax = self.y[self.max]
        valMax = np.ones(self.y.shape)
        valMin = np.zeros(self.y.shape)
        self.ny = (self.y - valMin) / (valMax-valMin)

    def bAct(self,b,c,A):
        b = np.asarray(b, order='c')
        c = np.asarray(c, order='c')
        A = np.asarray(A, order='c')
        return dot(dot(b,A),c.T).sum()

    def inverseCorrelation(self):
        K = eye(self.y.shape[0]) * self.obs_noise;
        for i,xa in enumerate(self.x):
            for j,xb in enumerate(self.x):
                K[i,j] += self.correlationFunction(xa,xb)

        return np.linalg.inv(K)
                

    def precomputeGPparams(self,Kinv):
        N = self.y.shape[0]
        uno = ones((N,))  #TODO: Generalize for other mean functions
        alpha = self.prior_alpha
        beta  = self.prior_beta
        delta = self.prior_delta
        uK = dot(uno,Kinv)
        eta = dot(uK,uno.T) + 1/delta
        yKy = self.bAct(self.ny,self.ny,Kinv)
        
        mu = dot(uK,self.ny) / eta;
        sigma2 = ( beta + yKy - mu*mu/eta ) / (alpha+N+2)

        return mu,sigma2,uK,eta


    def prediction(self,query,mu,sigma2,uK,eta,Kinv):
        N = self.y.shape[0]
        uno = ones((N,))  #TODO: Generalize for other mean functions
        r = array([])
        #TODO: try list comprehension
        for x in self.x:
            r = np.r_[r,self.correlationFunction(x,query)]

        rn = self.correlationFunction(query,query)

        rK = dot(r,Kinv)
        uKr = dot(uK,r)
        rKr = dot(rK,r)

        ymu = self.ny - mu

        y_pred  = mu + dot(rK,ymu.T)
        s_pred = np.sqrt(sigma2 + (rn - rKr + (1 - uKr)**2 / eta ))

        return y_pred,s_pred

    def plotResults(self):
        invR = self.inverseCorrelation();
        mu, sig2, uInvR, eta = self.precomputeGPparams(invR)
        xl = np.arange(0,1000) / 1000.0
        pl.figure()
        for x in xl:
            y,s = self.prediction(x,mu,sig2,uInvR,eta,invR)
            pl.plot(x,y,'k+')
            pl.plot(x,y+s,'r+')
            pl.plot(x,y-s,'r+')

        pl.plot(self.x,self.ny,'ko')

        pl.show()


    def run(self):
        for i in range(1,20):
            x = sp.rand()
            y = sp.rand()*50
            self.addNewPoint(x,y)

        self.plotResults()




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
    
