import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rc
plt.rc('figure', figsize=(6, 6))
from mpl_toolkits.mplot3d import Axes3D
import time

def telltime(seconds):
    h = np.floor(seconds / 3600)
    m = np.floor((seconds - 3600 * h) / 60)
    s = np.floor(seconds - 3600 * h - 60 * m)
    print(f"eta: {h} h {m} m {s} s")

# def AssignCells3d(cellcenters,positions,num):
#     particle2cell = np.ndarray(shape=(num, 3), dtype=int)
#     for i in range(len(positions)):
#         particle2cell[i,0] = np.argmin(np.abs(positions[i,0]-cellcenters))	# x coord of cell
#         particle2cell[i,1] = np.argmin(np.abs(positions[i,1]-cellcenters))	# y coord of cell
#         particle2cell[i,2] = np.argmin(np.abs(positions[i,2]-cellcenters))	# z coord of cell
#
#     return particle2cell
#
# def GenCell2Particle3d(particle2cell,num,ncells):
#     cell2particle = {}
#
#     for i in range(ncells):
#         for j in range(ncells):
#             for k in range(ncells):
#                 cell2particle[(i,j,k)] = []					# generates a list for each cell
#
#     for i in range(num):
#         cell2particle[tuple(particle2cell[i])].append(i) 	# assigns particles to lists
#
#     return cell2particle

np.random.seed(42)
N = 2
NT = 8*1400
L = 38.65
dt = 0.005

R = 1 													# sphere radius
diam = 2*R												# sphere diameter
dsq = diam*diam											# square sphere diameter

xx = np.zeros((N,NT),dtype=np.float64)			# particle positions x
yy = np.copy(xx)										# particle positions y
zz = np.copy(xx)										# particle positions z

vx = np.zeros((N,NT),dtype=np.float64)			# particle speeds x
vy = np.copy(vx)										# particle speeds y
vz = np.copy(vx)										# particle speeds z

# initialise speeds

#vx0 = np.random.rand(N)*2-1.0							# random initial speeds x (ranging -1 to 1)
#vy0 = np.random.rand(N)*2-1.0							# random initial speeds y (ranging -1 to 1)
#vz0 = np.random.rand(N)*2-1.0       					# random initial speeds z (ranging -1 to 1)

vx0 = np.array([1,1])
vy0 = np.array([1,0])
vz0 = np.array([1,0])

# initialise positions

x0 = np.array([0,-L/2])
y0 = np.array([-L/2,0])
z0 = np.array([-L/2,-L/2])

#x0 = np.random.rand(N)*L-0.5*L							# random initial position x (randing -L/2 to L/2)
#y0 = np.random.rand(N)*L-0.5*L							# random initial position y (randing -L/2 to L/2)
#z0 = np.random.rand(N)*L-0.5*L							# random initial position z (randing -L/2 to L/2)

xx[:,0] = x0											# initializes the x positions
yy[:,0] = y0											# initializes the y positions
zz[:,0] = z0											# initializes the z positions

vx[:,0] = vx0											# initializes the x speeds
vy[:,0] = vy0											# initializes the y speeds
vz[:,0] = vz0											# initializes the z speeds

# Define domain decomposition

ncells = int(np.floor(L/(2*R)))
cellsize = L/ncells
cellcenters = np.array([-0.5*L+(i+0.5)*cellsize for i in range(ncells)])

# Time-dependent adjacency matrix

am = np.zeros((N,N),dtype=np.float64)

step = 200
timestart = time.time()
for it in range(1,NT):

    if np.floor(it/step) > 0:
        print(f"completed update {it} out of {NT}")
        deltat = time.time() - timestart
        remainingtime = deltat/(200) * (NT - it)
        telltime(remainingtime)
        timestart = time.time()
        step += 200

    if np.abs( np.sum( np.abs(am)>1e-10 ) - 1.0*(N*N-N) ) < 1e-5:
        print("Adjacency matrix full, terminating early.")
        break # if adjacency matrix becomes full (apart from diagonal, kill simulation)

    # Update states (eventually, stop saving states, only save adjacency matrix)
    ###################################
    xx[:,it] = xx[:,it-1]+vx[:,it-1]*dt 	# Euler step for x
    yy[:,it] = yy[:,it-1]+vy[:,it-1]*dt		# Euler step for y
    zz[:,it] = zz[:,it-1]+vz[:,it-1]*dt		# Euler step for z

    vx[:,it] = vx[:,it-1]					# inertia x
    vy[:,it] = vy[:,it-1]					# inertia y
    vz[:,it] = vz[:,it-1]					# inertia z
    ###################################

    # implement reflecting BCs
    # checks if a position comp oversteps, if yes * corresponding v comp by -1
    ###################################
    # for i in range(N):
    #     if np.abs(xx[i,it]) > L/2:
    #         vx[i,it] *= -1
    #     if np.abs(yy[i,it]) > L/2:
    #         vy[i,it] *= -1
    #     if np.abs(zz[i,it]) > L/2:
    #         vz[i,it] *= -1
    ###################################

    # implement periodic BCs
    # checks if a position comp oversteps, if yes translates by L
    ###################################
    for i in range(N):
        if np.abs(xx[i, it]) > L / 2:
                xx[i, it] += -np.sign(xx[i, it])*L
        if np.abs(yy[i, it]) > L / 2:
                yy[i, it] += -np.sign(yy[i, it])*L
        if np.abs(zz[i, it]) > L / 2:
                zz[i, it] += -np.sign(zz[i, it])*L
    ###################################

    # check for "collisions" (all-to-all)
    ###################################
    # for i in range(N-1): # 0,...,N-2
    # 	for j in range(i+1,N): # i+1,...,N-1

    # 		sqdist = (xx[i,it]-xx[j,it])**2 + (yy[i,it]-yy[j,it])**2

    # 		# Should properly check that am is empty, and *only* record in that case
    # 		if sqdist < dsq and np.abs(am2[i,j])<1e-8:
    # 			am2[i,j] = it*dt
    # 			am2[j,i] = it*dt
    ###################################

    # check for "collisions" (domain decomposition)
    ###################################

    # matrix which takes particle (row) and returns cell index in x and y

    particle2cell = np.zeros(shape=(N, 3), dtype=np.int32)

    for i in range(N):
        particle2cell[i,0] = np.argmin(np.abs(xx[i,it] - cellcenters))
        particle2cell[i,1] = np.argmin(np.abs(yy[i,it] - cellcenters))
        particle2cell[i,2] = np.argmin(np.abs(zz[i,it] - cellcenters))

    # dictionary which takes (cellx, celly, cellz) and returns all particles in that cell

    cell2particle = {}
    for i in range(len(cellcenters)):
        for j in range(len(cellcenters)):
            for k in range(len(cellcenters)):
                cell2particle[(i,j,k)] = [] 		# probably faster to replace with: i + j*N "column/row-major indexing"
    for i in range(N):
        cell2particle[tuple(particle2cell[i])].append(i)

    for i in range(N):
        for iDx in [-1,0,1]:
            for iDy in [-1,0,1]:
                for iDz in [-1,0,1]:
                    cellx, celly, cellz = particle2cell[i]

                    cellx = cellx + iDx
                    celly = celly + iDy
                    cellz = cellz + iDz
                    if cellx >= 0 and cellx < ncells and celly >= 0 and celly < ncells and cellz >=0 and cellz < ncells:

                        # print(f"Particle {i}, cellx = {cellx}, celly = {celly}, cell2particle = {cell2particle[(cellx,celly)]}")
                        neighbors = cell2particle[(cellx,celly,cellz)]
                        for j in neighbors:
                            if j > i:
                                sqdist = (xx[i,it]-xx[j,it])**2 + (yy[i,it]-yy[j,it])**2 +(zz[i,it]-zz[j,it])**2

                                if sqdist < dsq and np.abs(am[i,j])<1e-8:
                                    am[i,j] = it*dt
                                    am[j,i] = it*dt
    ###################################

# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# ax.set_xlim(-L/2, L/2)
# ax.set_ylim(-L/2, L/2)
# ax.set_zlim(-L/2, L/2)
# ax.set_aspect('equal')
# for i in range(N):
#     plt.plot(xx[i].T,yy[i].T,zz[i].T)
# plt.show()
# print(np.max(np.abs(am-am.T)))
# print(am)

np.savetxt("collisionarray.csv", am, delimiter=",", fmt="%f")