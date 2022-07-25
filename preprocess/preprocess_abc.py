import pyvista as pv
import numpy as np
from scipy.spatial import Delaunay
import matplotlib.pyplot as plt
import sys

def compute_velocity(t, x, y, z):
    '''
        returns the velocity field of an ABC flow
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

def compute_velocity_vector(timeVector, points):
    velocity = []
    for i in range(len(timeVector)):
        for j in range(len(points)):
            velocity.append(compute_velocity(
                timeVector[i], points[j][0], points[j][1], points[j][2]))
    return velocity


# Generate space data
x = np.linspace(0, 1, 5)
xx, yy, zz = np.meshgrid(x, x, x)

# Set mesh points
points = np.column_stack((xx.ravel(order="F"),
                          yy.ravel(order="F"),
                          zz.ravel(order="F")))

# Create the point cloud mesh to triangulate from the coordinates
cloud = pv.PolyData(points)
mesh = cloud.delaunay_3d()
mesh = mesh.point_data_to_cell_data()
mesh.plot(show_edges=True)
n_cells = mesh.n_cells

tetrahedra = np.empty((n_cells, 4))
for i in range(n_cells):
    cell = mesh.extract_cells(i)
    # Find index of closest point in this mesh to the given point.
    tetrahedra[i][0] = mesh.find_closest_point(cell.points[0])
    tetrahedra[i][1] = mesh.find_closest_point(cell.points[1])
    tetrahedra[i][2] = mesh.find_closest_point(cell.points[2])
    tetrahedra[i][3] = mesh.find_closest_point(cell.points[3])

# Set initial time instants (velocity known in them)
time = np.linspace(0, 50, 100)
print(time)

# Set mesh velocity info
velocity = compute_velocity_vector(time, points)

# Store velocity, coords, time in files for C code
fc = open("abc_coords.txt", "w")
ff = open("abc_faces.txt", "w")
ft = open("abc_times.txt", "w")
fv = open("abc_velocity.txt", "w")

fc.write(str(len(points))+'\n')
ff.write(str(len(tetrahedra))+'\n')
ft.write(str(len(time))+'\n')

for index in range(len(points)):
    fc.write(str(points[index][0])+'\n')
    fc.write(str(points[index][1])+'\n')
    fc.write(str(points[index][2])+'\n')

for index in range(len(tetrahedra)):
    ff.write(str(tetrahedra[index][0])+'\n')
    ff.write(str(tetrahedra[index][1])+'\n')
    ff.write(str(tetrahedra[index][2])+'\n')
    ff.write(str(tetrahedra[index][3])+'\n')

for index in range(len(time)):
    ft.write(str(time[index])+'\n')

for index in range(len(velocity)):
    fv.write(str(velocity[index][0])+'\n')
    fv.write(str(velocity[index][1])+'\n')
    fv.write(str(velocity[index][2])+'\n')

# Close files
fc.close()
ff.close()
ft.close()
fv.close()
