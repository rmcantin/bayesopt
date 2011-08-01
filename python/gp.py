#!/usr/bin/env python

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

    def exponential(self,x1,x2):
        theta = self.params['theta']
        p = self.params['p']
        r = abs(x1-x2)
        dif = r**p / theta
        return np.exp(-dif.sum())

    def gaussian(self,x1,x2):
        theta = self.params['theta']
        r = x1-x2
        dif = r**2 / theta
        return np.exp(-dif.sum())
    
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

        self.obs_noise = 1e-4
        self.k = Kernel()
        self.normalize = False

    
    def correlationFunction(self,x1,x2):
        return self.k.mattern(x1,x2)

    def addNewPoint(self,newX,newY):
        self.x = np.r_[self.x,newX]
        self.y = np.r_[self.y,newY]

        minIndex = self.y.argmin()
        maxIndex = self.y.argmax()

        if self.normalize:
            normalize = False
            
            if minIndex != self.min:
                self.min = minIndex
                normalize = True

            if maxIndex != self.max:
                self.max = maxIndex
                normalize = True

            if normalize:
                self.normalizeResponseData()
        else:
            self.ny = np.r_[self.ny,newY]

    def normalizeResponseData(self):
        valMin = self.y[self.min]
        valMax = self.y[self.max]
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
                

    def fitGP(self):
        N = self.y.shape[0]
        uno = ones((N,))  #TODO: Generalize for other mean functions

        alpha = self.prior_alpha
        beta  = self.prior_beta
        delta = self.prior_delta
        
        Kinv = self.inverseCorrelation()
        uK = dot(uno,Kinv)
        eta = dot(uK,uno.T) + 1/delta

        yKy = self.bAct(self.ny,self.ny,Kinv)
        
        mu = dot(uK,self.ny) / eta;
        sigma2 = ( beta + yKy - mu*mu/eta ) / (alpha+N+2)

        return mu,sigma2,Kinv,uK,eta


    def prediction(self,query,mu,sigma2,Kinv,uK,eta):
        N = self.y.shape[0]
        uno = ones((N,))  #TODO: Generalize for other mean functions
        r = array([self.correlationFunction(x,query) for x in self.x])
        #rn = self.correlationFunction(query,query)
        rn = 1.0
        
        rK = dot(r,Kinv)
        uKr = dot(uK,r)
        rKr = dot(rK,r)

        ymu = self.ny - mu

        y_pred  = mu + dot(rK,ymu.T)
        s_pred = np.sqrt(sigma2 * (rn - rKr + (1 - uKr)**2 / eta ))

        return y_pred,s_pred

    def plotResults(self):
        mu,sig2,invR,uInvR,eta = self.fitGP()
        xl = np.arange(0,1000) / 1000.0
        pl.figure()
        for x in xl:
            y,s = self.prediction(x,mu,sig2,invR,uInvR,eta)

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

    def latinHypercubeSampling(self,nSamples,nDims):
        result = zeros(nSamples,nDims)
        f_nS = float(nSamples)
        
        for i in range(0,nDims):
            perms = np.random.permutation(nSamples)
            for j in range(0,nSamples):
                result[j,i] = (float(perms[j]) - np.random.rand()) / f_nS

        return result

    def uniformSampling(self,nSamples,nDims):
        return np.random.rand(nSamples,nDims)
        

    def evaluateQuery(self, query):
        if self.reachable(query):
            y,s = self.gp.prediction(query)
            if self.params_SKO_method == 'EI':
                return self.negEI(y,s)
            elif self.params_SKO_method == 'LCB':
                return self.lcb(y,s)
        
        return 0

    def lcb(self,y,s):
        a = 1.0
        return y - a*s

    def negEI(self,y,s):
        #Since data is normalized, ymin = 0
        yNorm = -y/s

        return y * norm.cdf(yNorm) - s * norm.pdf(yNorm)

    def negGEI(self,y,s,g):
        yN = -y/s
        pdfy = norm.pdf(yN)
        Tm2 = norm.cdf(yN)

        sumEI = (yN**g)*Tm2 - g*Tm1*yN**(g-1)

        for i in range(2,g):
            Tact = (i-1)*Tm2 - pdfy * yN**(i-1)
#            sumEI += -1**i * 

if __name__ == '__main__':
    myopt = GaussianProcess()
    myopt.run()
    
