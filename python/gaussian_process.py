#!/usr/bin/env python

import numpy as np
import pylab as pl
import scipy as sp
from numpy import eye,dot

def mattern(x1,x2,theta):
    # theta = sqrt(3) / l
    r = abs(x1 - x2)
    dif = r / theta
    acsum = dif.sum()
    acprod = np.prod(1+dif)

    return acprod*np.exp(-acsum)

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
        theta = 0.1
        return mattern(x1,x2,theta)
    
    def inverseCorrelation(self):
        obs_noise = 0.1
        K = eye(self.y.shape[0]) * obs_noise;
        for i,xa in enumerate(self.x):
            for j,xb in enumerate(self.x):
                K[i,j] += self.correlationFunction(xa,xb)

        return np.linalg.inv(K), np.linalg.det(K)

    def fitGP(self,x,y):
        self.x = x
        self.y = y

        N = y.shape[0]

        self.iK, detK = self.inverseCorrelation();
        self.lik = 0.5*(dot(dot(y,self.iK),y) - np.log(detK) - N*np.log(2*np.pi))

    def prediction(self,query):
        kv = np.array([self.correlationFunction(x,query) for x in self.x])
        kn = self.correlationFunction(query,query)

        kiK = dot(kv,self.iK)
        
        y_pred = dot(kiK,self.y.T)
        s_pred = kn - dot(kiK,kv)

        return y_pred,s_pred

    def plotResults(self):
        xl = np.arange(0,1000) / 1000.0
        pl.figure()
        for x in xl:
            y,s = self.prediction(x)

            pl.plot(x,y,'k+')
            pl.plot(x,y+s,'r+')
            pl.plot(x,y-s,'b+')

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
