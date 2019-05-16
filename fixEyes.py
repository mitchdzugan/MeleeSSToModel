import time
start = time. time()
import math
import os
import numpy as np
import cv2
import trimesh
import colorsys

# given a position (a, b, c) on a globe centered at the origin
# find the vector that represents the longitude line at that
# point in the direction of the north pole
def getLongRaw(a, b, c):
    z = 1 * math.pow((a*a) + (b*b), 1.0/2.0)
    z = z if c > 0 else -1 * z
    y = (-2*b*c*z) / ((2*a*a) + (2*b*b))
    inner = (c*c) - (y*y)
    x = math.sqrt(inner)
    x = x if a < 0 else -1 * x
    return [x, y, z]

# magniude of a vector
def mag(x, y, z):
    return math.sqrt(x*x + y*y + z*z)

# return a vector with same direction with magniude 1
def norm(x, y, z):
    mag_ = mag(x, y, z)
    return [x/mag_, y/mag_, z/mag_]

# cross product
def cross(Ax, Ay, Az, Bx, By, Bz):
    x = (Ay*Bz) - (Az*By)
    y = (Az*Bx) - (Ax*Bz)
    z = (Ax*By) - (Ay*Bx)
    return norm(x, y, z)

# dot product
def dot(A, B):
    return A[0]*B[0] + A[1]*B[1] + A[2]*B[2]

# gets vectors representing the x and y axis
# for a picture taken from point (a, b, c)
# directed at the origin
def getAxis(a, b, c):
    [Yx, Yy, Yz] = norm(*getLongRaw(a, b, c))
    [Xx, Xy, Xz] = cross(Yx, Yy, Yz, a, b, c)
    return [[Xx, Xy, Xz], [Yx, Yy, Yz]]


def getDimFromFilename(fname):
    res = {}
    res['fname'] = fname
    fname = fname.replace('.png', '')
    pieces = fname.split('_')
    print(pieces)
    res['tx'] = float(pieces[0])
    res['ty'] = float(pieces[1])
    res['tz'] = float(pieces[2])
    res['cx'] = float(pieces[3])
    res['cy'] = float(pieces[4])
    res['cz'] = float(pieces[5])
    return res

STEP = 1.0/30.0

XRANGE = int(14.5 / STEP)
YRANGE = int(17.0 / STEP)
ZRANGE = int(11.0 / STEP)

XSTART = -11.7
YSTART = 51.0
ZSTART = -6.0

matrix = np.ones((XRANGE, YRANGE, ZRANGE), dtype=np.bool)

for x_step in range(XRANGE):
    for y_step in range(YRANGE):
        for z_step in range(ZRANGE):
            matrix[x_step, y_step, z_step] = True

matrix = np.load('falco.all.npy')

def getPixel(img, x, y):
    x = int(x)
    y = int(y)
    if (x < 0 or y < 0):
        return [[0, 255, 0], False]
    if (x >= img.shape[1] or y >= img.shape[0]):
        return [[0, 255, 0], False]
    return [img[y, x], True]

# find the x, y position in the screenshot
# where this pos [x, y, z] would appear
def getPicPos(img, dims, normView, axis, pos):
    isUnder = axis[1][2] < 0
    centerX = img.shape[1] / 2
    centerY = img.shape[0] / 2
    [x, y, z] = pos

    fromCam = [dims['cx']-x, dims['cz']-z, dims['cy']-y]
    dist = dot(fromCam, normView)
    offsetX = dot(fromCam, axis[0])
    offsetY = dot(fromCam, axis[1])

    upp = 0.00102454 * dist
    ppu = 1/(upp)

    offsetX = (offsetX * ppu * img.shape[1] / 636.0)
    offsetY = (offsetY * ppu * img.shape[0] / 524.0)
    offsetX = offsetX if not isUnder else -1 * offsetX
    offsetY = offsetY if not isUnder else -1 * offsetY
    picX = centerX + offsetX
    picY = centerY + offsetY

    return [picX, picY, dist]


def makeKey(x, y):
    return `x` + '_' + `y`
def processEyes(dname, fname, pxmin, pxmax, pymin, pymax):
    img = cv2.imread(dname + '/' + fname)

    totalDist = 0
    dims = getDimFromFilename(fname)
    normView = norm(dims['cx']-dims['tx'], dims['cz']-dims['tz'], dims['cy']-dims['ty'])
    [x_axis, y_axis] = getAxis(*normView)


    eyeDists = {}
    print(x_axis)
    print(y_axis)
    for x_step in range(XRANGE):
        print([x_step, XRANGE])
        for y_step in range(YRANGE):
            for z_step in range(ZRANGE):
                curr = matrix[x_step, y_step, z_step]
                if curr:
                    x = XSTART + x_step * STEP
                    y = YSTART + y_step * STEP
                    z = ZSTART + z_step * STEP

                    picPos = getPicPos(
                        img,
                        dims,
                        normView,
                        [x_axis, y_axis],
                        [x, y, z]
                    )

                    [picX, picY, dist] = picPos
                    picX = int(round(picX))
                    picY = int(round(picY))
                    if picX >= pxmin and picX <= pxmax and picY >= pymin and picY <= pymax:
                        [b, g, r] = img[picY, picX]
                        [h, s, v] = colorsys.rgb_to_hsv(r / 255.0, g / 255.0, b / 255.0)
                        if (h > 0.3888 and h < 0.6111) or s < 0.15 or v < 0.15:
                            currKey = makeKey(picX, picY)
                            if eyeDists.has_key(currKey):
                                eyeDists[currKey] = min(dist, eyeDists[currKey])
                            else:
                                eyeDists[currKey] = dist
    for x_step in range(XRANGE):
        print([x_step, XRANGE])
        for y_step in range(YRANGE):
            for z_step in range(ZRANGE):
                curr = matrix[x_step, y_step, z_step]
                if curr:
                    x = XSTART + x_step * STEP
                    y = YSTART + y_step * STEP
                    z = ZSTART + z_step * STEP

                    picPos = getPicPos(
                        img,
                        dims,
                        normView,
                        [x_axis, y_axis],
                        [x, y, z]
                    )

                    [picX, picY, dist] = picPos
                    picX = int(round(picX))
                    picY = int(round(picY))
                    currKey = makeKey(picX, picY)
                    if eyeDists.has_key(currKey):
                        if eyeDists[currKey] + (5.0 * STEP) > dist:
                            matrix[x_step, y_step, z_step] = False
    return

dname = "/mnt/c/Users/Mitch/Pictures/FALCO_EYES_RUL"
fname1 = "-6.159353256225586_59.242515563964844_-5.290919303894043_11.40567398071289_59.121578216552734_37.49576950073242.png"
fname2 = "-6.159353256225586_59.242515563964844_-5.290919303894043_19.295011520385742_82.3686294555664_-36.21831512451172.png"
processEyes(dname, fname1, 533, 604, 273, 319)
processEyes(dname, fname2, 361, 407, 233, 289)

np.save('falco.all.eyes.npy', matrix)
mesh = trimesh.voxel.matrix_to_marching_cubes(
                matrix=matrix,
                pitch=1.0,
                origin=np.zeros(3))

mesh.export('falco.all.eyes.obj')
end = time. time()
print(end - start)
