#!/usr/bin/env python

import numpy as np
import pylab as pl
import scipy as sp
from numpy import eye,dot
from scipy import spatial

def mattern(x1,x2,params):

    #TODO: Use lambda functions to deal with different mattern
    # kernels d = params["d"];
    th = params["theta"]
    r = np.atleast_2d(np.r_[x1,x2])
    mh = np.sqrt(spatial.distance.pdist(r.T,'mahalanobis',VI=eye(r.shape[0])*th))
    #dif = r / theta
    k = (1+mh)*np.exp(-mh)
    dkdt = mh*mh*np.exp(-mh) 
    return k

def exponential(x1,x2,theta,p):
    r = abs(x1-x2)
    dif = r**p / theta
    return np.exp(-dif.sum())

def gaussian(x1,x2,theta):
    r = x1-x2
    dif = r**2 / theta
    return np.exp(-dif.sum())

class BasicGaussianProcess:
    
    def correlationFunction(self,x1,x2):
        params = {"theta" : 0.1}
        return mattern(x1,x2,params)
    
    def inverseCorrelation(self):
        obs_noise = 0.01
        n = self.y.shape[0]
        K = eye(n) * obs_noise
        dK = np.zeros((n,n))
        for i,xa in enumerate(self.x):
            for j,xb in enumerate(self.x):
                K[i,j], dK[i,j] += self.correlationFunction(xa,xb)

        return np.linalg.inv(K), dK, np.linalg.det(K)

    def fitGP(self,x,y):
        self.x = x
        self.y = y

        N = y.shape[0]

        self.iK, self.dK, detK = self.inverseCorrelation();
        av = dot(self.iK,self.y)
        
        lik = 0.5*(dot(dot(y,self.iK),y) - np.log(detK) - N*np.log(2*np.pi))
        grad = 0.5*numpy.trace(dot(dot(av,av.T) - self.iK,self.dK))

        return lik, grad

    def prediction(self,query):
        kv = np.array([self.correlationFunction(x,query) for x in self.x])
        kn = self.correlationFunction(query,query)

        kiK = dot(kv.T,self.iK)
        
        y_pred = np.atleast_2d(dot(kiK,self.y.T))
        s_pred = np.sqrt(kn - dot(kiK,kv))

        return y_pred,s_pred

    def plotResults(self):
        xl = np.arange(0,1000) / 1000.0
        yl = np.atleast_2d(np.array([])).T
        sl = np.atleast_2d(np.array([])).T
        pl.figure()
        for x in xl:
            y,s = self.prediction(x)
            yl = np.r_[yl,y]
            sl = np.r_[sl,s]
            

        pl.plot(xl,yl,'k-')
        pl.plot(xl,yl+3*sl,'r-')
        pl.plot(xl,yl-3*sl,'r-')

        pl.plot(self.x,self.y,'ko')

        pl.show()

    def generateSampleStepFunction(self):
        x = sp.rand()
        y = sp.randn()*0.1
        
        if x > 0.5:
            y += 10
        return x,y

    def run(self):
        xl = np.array([])
        yl = np.array([])
        for i in range(1,100):
            x,y = self.generateSampleStepFunction()
            xl = np.r_[xl,x]
            yl = np.r_[yl,y]
            
        self.fitGP(xl,yl)
        self.plotResults()

if __name__ == '__main__':
    myopt = BasicGaussianProcess()
    myopt.run()
