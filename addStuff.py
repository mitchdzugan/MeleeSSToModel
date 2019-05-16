import time
start = time. time()
import math
import os
import numpy as np
import cv2
import trimesh

STEP = 1.0/30.0

XRANGE = int(14.5 / STEP)
YRANGE = int(17.0 / STEP)
ZRANGE = int(11.0 / STEP)

XSTART = -6.5
YSTART = 51.0
ZSTART = -6.0

matrix = np.load('falco.all.base.17.npy')

ybegin = int(YRANGE * 0.59)
xbegin = int(XRANGE * 0.45)
xend = int(XRANGE * 0.25)
width = abs(xbegin - xend)
leahboo = cv2.imread('./leahboo.png')

scale = 1.0 * width / leahboo.shape[1]
height = scale * leahboo.shape[0]
print(leahboo.shape)
print(width)
print(height)

# add leahboo to jacket
for xstep in range(width):
    for ystep in range(int(height)):
        ximg = min(int(xstep / scale), leahboo.shape[1] - 1)
        yimg = min(int(ystep / scale), leahboo.shape[0] - 1)
        pix = leahboo[yimg, ximg]
        v = pix[0]
        if (v > 237):
            continue
        vv = v / (238.0 / 3)
        vvv = 3 + 3 - vv
        zmin = ZRANGE / 2
        print({'vvv': vvv, 'vv': vv, 'v': v})
        while True:
            if not matrix[xbegin - xstep, ybegin - ystep, zmin]:
                break
            zmin -= 1
        for zoffset in range(int(vvv)):
            matrix[xbegin - xstep, ybegin - ystep, zmin - zoffset] = True


# add base platform
x_mid = XRANGE / 2.0
z_mid = ZRANGE / 2.0
for x_step in range(XRANGE):
    print([x_step, XRANGE])
    for y_step in range(YRANGE):
        for z_step in range(ZRANGE):
            if y_step < 0.04 * YRANGE:
                if z_step > 0.20 * ZRANGE:
                    if x_step < 0.99 * XRANGE and x_step > 0.01 * XRANGE:
                        matrix[x_step, y_step, z_step] = True

np.save('falco.all.npy', matrix)
mesh = trimesh.voxel.matrix_to_marching_cubes(
                matrix=matrix,
                pitch=1.0,
                origin=np.zeros(3))

mesh.export('falco.all.obj')
end = time. time()
print(end - start)
