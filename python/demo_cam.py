#!/usr/bin/env python
# -------------------------------------------------------------------------
#    This file is part of BayesOpt, an efficient C++ library for 
#    Bayesian optimization.
#
#    Copyright (C) 2011-2015 Ruben Martinez-Cantin <rmcantin@unizar.es>
# 
#    BayesOpt is free software: you can redistribute it and/or modify it 
#    under the terms of the GNU Affero General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    BayesOpt is distributed in the hope that it will be useful, but 
#    WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Affero General Public License for more details.
#
#    You should have received a copy of the GNU Affero General Public License
#    along with BayesOpt.  If not, see <http://www.gnu.org/licenses/>.
# ------------------------------------------------------------------------

from SimpleCV import Camera
import numpy as np
import bayesopt
from time import sleep

# Initialize the camera
cam = Camera()
cost = np.zeros(256)

#Load images
img = cam.getImage().scale(400,400)
img2 = img.binarize()


def costImage(i):
    # Make image black and white
    img1 = img.binarize(int(i))
    mat = img1.getNumpy()
    countW = np.count_nonzero(mat);
    countB = mat.size-countW
    return ((countW-countB)/float(mat.size))**2

params = {} #bayesopt.initialize_params()
params['n_iterations'] = 15
params['n_init_samples'] = 5

valid_values = np.transpose(np.array(range(256), dtype=float, ndmin=2))
mvalue, x_out, error = bayesopt.optimize_discrete(costImage,
                                                  valid_values, params)

x_out = int(x_out)
print x_out
img1 = img.binarize(x_out)

img1 = img.sideBySide(img1).sideBySide(img2)
img1.drawText("Threshold: "+str(x_out))
img1.show()

foo = raw_input('Press any key')
