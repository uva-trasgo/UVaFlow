import pyvista as pv
import numpy as np
import sys
import math

def compute_velocity(t, x, y):
    A = 0.1
    omega = 2*np.pi/10
    epsilon = 0.25
    a_t = epsilon*np.sin(omega*t)
    b_t = 1 - 2*epsilon*np.sin(omega*t)
    f = a_t*x**2+b_t*x
    dfdx = 2*a_t*x+b_t
    u = -A*np.pi*np.sin(np.pi*f)*np.cos(np.pi*y)
    v = np.pi * A*np.cos(np.pi*f)*np.sin(np.pi*y)*dfdx
    return [u, v]

def compute_velocity_vector(timeVector, points):
    velocity = []
    for i in range(len(timeVector)):
        for j in range(len(points)):
            velocity.append(compute_velocity(timeVector[i], points[j][0], points[j][1]))
    return velocity

# Load Dolphin data
x = np.linspace(0, 2, 100)
y = np.linspace(0, 1, 50)
xx, yy, zz = np.meshgrid(x, y, [0])

points = np.column_stack((xx.ravel(order="F"),
                          yy.ravel(order="F"),
                          zz.ravel(order="F")))

cloud = pv.PolyData(points)
mesh = cloud.delaunay_2d()

time = np.linspace(0, 50, 501)
print(time)

for i, t_ in enumerate(time):
    t_=round(t_,1)
    if(t_-round(t_,0)==0):
        t_=int(t_)
    U, V = compute_velocity(mesh.points[:, 0], mesh.points[:, 1], t_)
    mesh['U'] = U
    mesh['V'] = V
    #filename = './c7_triangles_vtk/bottom_plane/'+str(t_)+'/zplane'
    #mesh.save(filename+'.vtp')

assert mesh.is_all_triangles()

# Extract points and faces
# Note: faces reshaped taking into account all are triangles (i.e. all have 3 vertices)
points = mesh.points
faces = mesh.faces.reshape((-1,4))[:, 1:4]

# Set mesh velocity info
velocity = compute_velocity_vector(time, points)

# Store velocity, coords, time in files for C code
fc = open("doublegire_coords.txt", "w")
ff = open("doublegire_faces.txt", "w")
ft = open("doublegire_times.txt", "w")
fv = open("doublegire_velocity.txt", "w")

fc.write(str(len(points))+'\n')
ff.write(str(len(faces))+'\n')
ft.write(str(len(time))+'\n')

for index in range(len(points)):
    fc.write(str(points[index][0])+'\n')
    fc.write(str(points[index][1])+'\n')

for index in range(len(faces)):
    ff.write(str(faces[index][0])+'\n')
    ff.write(str(faces[index][1])+'\n')
    ff.write(str(faces[index][2])+'\n')

for index in range(len(time)):
    ft.write(str(time[index])+'\n')

for index in range(len(velocity)):
    fv.write(str(velocity[index][0])+'\n')
    fv.write(str(velocity[index][1])+'\n')

# Close files
fc.close()
ff.close()
ft.close()
fv.close()
