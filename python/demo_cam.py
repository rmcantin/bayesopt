from SimpleCV import Camera
import numpy as np
# Initialize the camera
cam = Camera()
cost = np.zeros(256)

#Load images
img = cam.getImage().scale(200,200)
img2 = img.binarize()

# Loop to continuously get images
for i in range(256):
    # Get Image from camera

    # Make image black and white
    img1 = img.binarize(i)

    mat = img1.getNumpy()
    countW = np.count_nonzero(mat);
    countB = mat.size-countW
    cost[i] = (countW-countB)**2
    print countW, countB, cost[i]

    # Show the image
    img1 = img.sideBySide(img1).sideBySide(img2)
    img1.drawText("Threshold: "+str(i))
    img1.show()

print np.argmin(cost)
