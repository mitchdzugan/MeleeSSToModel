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

def isGreen(img, x, y):
    xl = x - 0.5
    xr = x + 0.5
    yt = y - 0.5
    yb = y + 0.5

    [pbl, onScreen1] = getPixel(img, xl, yb)
    [pbr, onScreen2] = getPixel(img, xr, yb)
    [ptl, onScreen3] = getPixel(img, xl, yt)
    [ptr, onScreen4] = getPixel(img, xr, yt)

    offScreen = not onScreen1 or not onScreen2 or not onScreen3 or not onScreen4

    wbl = (math.ceil(xl) - xl) * (math.ceil(yb) - yb)
    wbr = (xr - math.floor(xr)) * (math.ceil(yb) - yb)
    wtl = (math.ceil(xl) - xl) * (yt - math.floor(yt))
    wtr = (xr - math.floor(xr)) * (yt - math.floor(yt))

    point = [
        wbl*pbl[0] + wbr*pbr[0] + wtl*ptl[0] + wtr*ptr[0],
        wbl*pbl[1] + wbr*pbr[1] + wtl*ptl[1] + wtr*ptr[1],
        wbl*pbl[2] + wbr*pbr[2] + wtl*ptl[2] + wtr*ptr[2],
    ]

    return [point[0] < 5 and point[1] > 100 and point[2] < 5, offScreen]

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

def processFile(dname, fname, skipOffScreen):
    img = cv2.imread(dname + '/' + fname)

    totalDist = 0
    dims = getDimFromFilename(fname)
    normView = norm(dims['cx']-dims['tx'], dims['cz']-dims['tz'], dims['cy']-dims['ty'])
    [x_axis, y_axis] = getAxis(*normView)

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

                    [isCurrGreen, isOffScreen] = isGreen(img, picX, picY)
                    if (skipOffScreen and (isOffScreen or dist < 0)):
                        pass
                    else:
                        matrix[x_step, y_step, z_step] = not isCurrGreen
    return


def doRun(dname, skipOffScreen):
    pics = os.listdir(dname)

    i_xd = 0
    for pic in pics:
        print(['PIC', dname, pic, i_xd, len(pics)])
        i_xd += 1
        processFile(dname, pic, skipOffScreen)

# doRun('/mnt/c/Users/Mitch/Pictures/FALCO_ENDGAME/CORE', False)
# doRun('/mnt/c/Users/Mitch/Pictures/FALCO_ENDGAME/PARTIAL_15', True)
# doRun('/mnt/c/Users/Mitch/Pictures/FALCO_ENDGAME/UNDER', False)


np.save('falco.all.npy', matrix)
mesh = trimesh.voxel.matrix_to_marching_cubes(
                matrix=matrix,
                pitch=1.0,
                origin=np.zeros(3))

mesh.export('falco.all.obj')
end = time. time()
print(end - start)

xmin = None
xmax = None
for x_step in range(XRANGE):
    for y_step in range(YRANGE):
        for z_step in range(ZRANGE):
            curr = matrix[x_step, y_step, z_step]
            if curr:
                x = XSTART + x_step * STEP
                if xmin == None or xmin > x:
                    xmin = x
                if xmax == None or xmax < x:
                    xmax = x

print([xmin, xmax])
