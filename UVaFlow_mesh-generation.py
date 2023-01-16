#! /usr/bin/env python3

import pyvista as pv
import numpy as np
import time
import sys as sys
from joblib import Parallel, delayed
import multiprocessing
from multiprocessing import shared_memory
from scipy.spatial import cKDTree as KDTree

import argparse
from ctypes import *

from scipy.interpolate import RegularGridInterpolator
from scipy.integrate import solve_ivp

def compute_velocity_vector_2D_DGyre(timeVector, points):
    velocity = []
    for i in range(len(timeVector)):
        for j in range(len(points)):
            velocity.append(compute_velocity_2D_DGyre(timeVector[i], points[j][0], points[j][1]))
    return velocity

def compute_velocity_2D_DGyre(t, x, y):
    '''
        returns the velocity field of a Double Gyre flow (2D)
    '''
    A = 0.1
    omega = 2*np.pi/10
    epsilon = 0.25
    a_t = epsilon*np.sin(omega*t)
    b_t = 1 - 2*epsilon*np.sin(omega*t)
    f = a_t*x**2+b_t*x
    dfdx = 2*a_t*x+b_t
    u = -A*np.pi*np.sin(np.pi*f)*np.cos(np.pi*y)
    v = np.pi * A*np.cos(np.pi*f)*np.sin(np.pi*y)*dfdx
    return u, v

def compute_velocity_3D_ABC(t, x, y, z):
    '''
        returns the velocity field of an ABC flow (3D)
    '''
    A = np.sqrt(3)
    B = np.sqrt(2)
    C = 1
    omega = 2*np.pi/10
    epsilon = 0.1
    A_t = A+epsilon*np.cos(omega*t)
    u = A_t*np.sin(z) + C*np.cos(y)
    v = B*np.sin(x) + A_t*np.cos(z)
    w = C*np.sin(y) + B*np.cos(x)

    return u, v, w

def compute_velocity_vector_3D_ABC(timeVector, points):
    velocity = []
    for i in range(len(timeVector)):
        for j in range(len(points)):
            velocity.append(compute_velocity_3D_ABC(
                timeVector[i], points[j][0], points[j][1], points[j][2]))
    return velocity

def main_2D_DGyre (args):
    
    t_start = time.time()

    # Parse args
    nDim = int(args.nDim[0])
    coords_file = args.coords_file
    faces_file = args.faces_file
    times_file = args.times_file
    vel_file = args.vel_file

    print("x_steps_axis -> %d" % int(args.x_steps_axis))
    print("y_steps_axis -> %d" % int(args.y_steps_axis))

    x_ = np.linspace(0, 2, int(args.x_steps_axis))
    y_ = np.linspace(0, 1, int(args.y_steps_axis))
    t_ = np.linspace(float(args.t_zero), float(args.t_end), int(args.nt_knwon))
    xx, yy, zz = np.meshgrid(x_, y_, [0])

    points = np.column_stack((xx.ravel(order="F"),
                            yy.ravel(order="F"),
                            zz.ravel(order="F")))
    print("Number of points: %d" % len(points), flush=True)

    cloud = pv.PolyData(points)
    print("Constructing 2D Delaunay triangulation of the mesh...", flush=True)
    mesh = cloud.delaunay_2d()

    xg, yg, tg = np.meshgrid(x_, y_, t_, indexing='ij')
    _v = np.array(compute_velocity_2D_DGyre(tg, xg, yg))
    vx = _v[0]
    vy = _v[1]

    gridpoints = np.array(mesh.points)
    vel = compute_velocity_vector_2D_DGyre(t_, gridpoints)

    gridpoints = np.delete(gridpoints, 2, 1)

    # Store velocity, coords, time in files for the UVaFlow C code
    t_preprocess = time.time()
    print("Python data generated: "+str(t_preprocess-t_start), flush=True)

    fc = open(coords_file, "w")
    ff = open(faces_file, "w")
    ft = open(times_file, "w")
    fv = open(vel_file, "w")

    fc.write(str(len(points))+'\n')
    for index in range(len(points)):
        fc.write(str(points[index][0])+'\n')
        fc.write(str(points[index][1])+'\n')

    ft.write(str(len(t_))+'\n')
    for index in range(len(t_)):
       ft.write(str(t_[index])+'\n')
       for jindex in range(len(points)):
           fv.write(str(vel[index*len(points)+jindex][0])+'\n')
           fv.write(str(vel[index*len(points)+jindex][1])+'\n')

    faces = mesh.faces.reshape((-1,4))[:, 1:4]
    print("Number of faces: %d" % len(faces) , flush=True )
    ff.write(str(len(faces))+'\n')    
    for index in range(len(faces)):
        ff.write(str(faces[index][0])+'\n')
        ff.write(str(faces[index][1])+'\n')
        ff.write(str(faces[index][2])+'\n')

    # Close files
    fc.close()
    ff.close()
    ft.close()
    fv.close()

    t_store = time.time()
    print("Python data saved in files: "+str(t_store-t_preprocess), flush=True)

    t_end = time.time()
    print("Total time elapsed: "+str(t_end-t_start), flush=True)

def main_3D_ABC (args):

    t_start = time.time()

    # Parse args
    nDim = int(args.nDim[0])
    coords_file = args.coords_file
    faces_file = args.faces_file
    times_file = args.times_file
    vel_file = args.vel_file

    print("x_steps_axis -> %d" % int(args.x_steps_axis))
    print("y_steps_axis -> %d" % int(args.y_steps_axis))
    print("z_steps_axis -> %d" % int(args.z_steps_axis))

    x_ = np.linspace(0, 1, int(args.x_steps_axis))
    y_ = np.linspace(0, 1, int(args.y_steps_axis))
    z_ = np.linspace(0, 1, int(args.z_steps_axis))
    t_ = np.linspace(float(args.t_zero), float(args.t_end), int(args.nt_knwon))
    xx, yy, zz = np.meshgrid(x_, y_, z_)

    points = np.column_stack((xx.ravel(order="F"),
                            yy.ravel(order="F"),
                            zz.ravel(order="F")))

    print("Number of points: %d" % len(points), flush=True)
    cloud = pv.PolyData(points)
    print("Constructing 3D Delaunay triangulation of the mesh...", flush=True)
    mesh = cloud.delaunay_3d()
    mesh = mesh.point_data_to_cell_data()
    n_cells = mesh.n_cells
    
    data = np.empty((n_cells, 4))
    shm = shared_memory.SharedMemory(name='shared_tetrahedra3', create=True, size=data.nbytes)
    tetrahedra = np.ndarray(data.shape, dtype=str(data.dtype), buffer=shm.buf)
    
    print("Generating KDTree...", flush=True)
    # KDTree
    kdtree_mesh_points = KDTree(mesh.points)

    # Ensure mesh and data exists
    def define_tetrahedra(p):
        ini = p[0]
        fin = p [1]
        t_local = []
        for i in range(fin - ini):
            cell = mesh.extract_cells(ini + i)
            # Get shared memory
            distances, indices = kdtree_mesh_points.query(cell.points)
            t_local.append(indices)

        existing_shm = shared_memory.SharedMemory(name='shared_tetrahedra3')
        t = np.ndarray(data.shape, dtype=str(data.dtype), buffer=existing_shm.buf)
        #Find index of closest point in this mesh to the given point.
        for i in range(fin - ini):
            i_ = ini+i
            t[i_][0] = t_local[i][0]
            t[i_][1] = t_local[i][1]
            t[i_][2] = t_local[i][2]
            t[i_][3] = t_local[i][3]
        existing_shm.close()

    print("Generating tetrahedra...", flush=True)

    # Create tetrahedra in parallel
    num_cores = int(args.num_cores)

    step = n_cells // num_cores
    start = 0
    stop = step
    parts = []
    for i in range(num_cores):
        parts.append( (start, stop) )
        start += step
        stop += step
        #last core
        if i+1 == (num_cores-1):
            stop += n_cells % num_cores

    results = Parallel(n_jobs=num_cores)(delayed(define_tetrahedra)(p) for p in parts)

    xg, yg, zg, tg = np.meshgrid(x_, y_, z_, t_, indexing='ij')
    _v = np.array(compute_velocity_3D_ABC(tg, xg, yg, zg) )
    vx = _v[0]
    vy = _v[1]
    vz = _v[2]

    gridpoints = np.array(mesh.points)
    vel = compute_velocity_vector_3D_ABC(t_, gridpoints)

    t_preprocess = time.time()
    print("Python data generated: "+str(t_preprocess-t_start), flush=True)

    fc = open(coords_file, "w")
    ff = open(faces_file, "w")
    ft = open(times_file, "w")
    fv = open(vel_file, "w")

    fc.write(str(len(points))+'\n')
    for index in range(len(points)):
        fc.write( str(points[index][0])+'\n'+ str(points[index][1])+'\n' + str(points[index][2])+'\n')

    ft.write(str(len(t_))+'\n')
    for index in range(len(t_)):
       ft.write(str(t_[index])+'\n')
       for jindex in range(len(points)):
           fv.write(str(vel[index*len(points)+jindex][0])+'\n'+str(vel[index*len(points)+jindex][1])+'\n'+str(vel[index*len(points)+jindex][2])+'\n')

    ff.write(str(len(tetrahedra))+'\n')
    for index in range(len(tetrahedra)):
        ff.write(str(tetrahedra[index][0])+'\n' + str(tetrahedra[index][1])+'\n' + str(tetrahedra[index][2])+'\n' + str(tetrahedra[index][3])+'\n')

    # Shared memory end
    shm.close()
    shm.unlink()

    # Close files
    fc.close()
    ff.close()
    ft.close()
    fv.close()

    t_store = time.time()
    print("Python data saved in files: "+str(t_store-t_preprocess), flush=True)
    t_end = time.time()
    print("Total time elapsed: "+str(t_end-t_start), flush=True)

if __name__ == '__main__':

    #COMPUTE_FLOWMAP_C_FUNCTIONS = CDLL(SHARED_LIBRARY)

    NDIM_CHOICES = ["2D_DGyre","3D_ABC"]
    NDIM_HELP = "Dimensions of the space. Choices: %s" % (", ".join(NDIM_CHOICES) )
    COMPUTE_FLOWMAP_POLICIES = ["SEQUENTIAL", "OMP_STATIC", "OMP_DYNAMIC", "OMP_GUIDED"]
    COMPUTE_FLOWMAP_POLICIES_HELP = "Scheduling policy. Choices: %s" % (", ".join(COMPUTE_FLOWMAP_POLICIES))

    parser = argparse.ArgumentParser()
    parser.add_argument(metavar="<nDim>", dest="nDim", help=NDIM_HELP, choices=NDIM_CHOICES, default="3D_ABC")
    parser.add_argument(metavar="<nt_knwon>", dest="nt_knwon", help="Number of time instants for which the velocity is knwon")
    parser.add_argument(metavar="<t_zero>", dest="t_zero", help="First time instant for which the velocity is knwon")
    parser.add_argument(metavar="<t_end>", dest="t_end", help="Last time instant for which the velocity is knwon")
    parser.add_argument(metavar="<coords_file>", dest="coords_file", help="File path where mesh coordinates will be stored", default="coords.txt")
    parser.add_argument(metavar="<faces_file>", dest="faces_file", help="File path where mesh faces will be stored", default="faces.txt")
    parser.add_argument(metavar="<times_file>", dest="times_file", help="File path where time data will be stored", default="times.txt")
    parser.add_argument(metavar="<vel_file>", dest="vel_file", help="File path where velocity data will be stored", default="velocity.txt")
    parser.add_argument(metavar="<num_cores>", dest="num_cores", help="Threads to use", type=int)
    parser.add_argument(metavar="<x_steps_axis>", dest="x_steps_axis", help="Steps in X axis (for linspace)",  default=-1, type=int )
    parser.add_argument(metavar="<y_steps_axis>", dest="y_steps_axis", help="Steps in Y axis (for linspace)", default=-1, type=int )
    parser.add_argument(metavar="<z_steps_axis>", dest="z_steps_axis", help="Steps in Z axis (for linspace)", default=-1, nargs='?', type=int)

    args = parser.parse_args()
    
    # Ensure z_steps_axis is defined for 3D
    if args.x_steps_axis <= 0 or args.y_steps_axis <= 0:
        print("ERROR: <x_steps_axis> and <y_steps_axis> must be positive integers")
        sys.exit()
    if args.nDim == "3D" and args.z_steps_axis == -1:
        print("ERROR: you must define <z_steps_axis> for 3D")
        parser.print_help()
        sys.exit()

    eval("main_"+args.nDim)(args)
