import time
start = time. time()
import math
import os
import numpy as np
import cv2
import trimesh

import sys
sys.setrecursionlimit(10**9)

# removes any lose parts of the mesh that are
# not connected to the main body
matrix = np.load('falco.all.eyes.npy')
mesh = trimesh.voxel.matrix_to_marching_cubes(
                matrix=matrix,
                pitch=1.0,
                origin=np.zeros(3))
meshes = mesh.split()
meshes[0].export('falco.complete.eyes.obj')
end = time. time()
print(end - start)
