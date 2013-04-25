from SimpleCV import Camera
import numpy as np
import bayesopt
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

params = bayesopt.initialize_params()
params['n_iterations'] = 15
params['n_init_samples'] = 10
params['sigma_s'] = 1
params['crit_name'] = "cEI"
params['kernel_name'] = "kMaternISO3"
params['mean_name'] = "mOne"

valid_values = np.transpose(np.array(range(256), dtype=float, ndmin=2))
mvalue, x_out, error = bayesopt.optimize_discrete(costImage,
                                                  valid_values, params)

x_out = int(x_out)
print x_out
img1 = img.binarize(x_out)

# Loop to continuously get images
for i in range(256):
    cost[i] = costImage(i)

minid = np.argmin(cost)

print minid, cost[minid]
img3 = img.binarize(np.argmin(cost))

img1 = img.sideBySide(img1).sideBySide(img2)
img1.drawText("Threshold: "+str(x_out))
img1.show()

foo = raw_input('Press any key')
