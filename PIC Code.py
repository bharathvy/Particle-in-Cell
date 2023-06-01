from numpy import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import time

# Defining simulation parameters

Np = 50000  # number of particles
Lx = 0.1    # length of simulation box in x direction
Ly = 0.1    # length of simulation box in y direction
Lz = 0.1    # length of simulation box in z direction
dx = 0.01     # grid spacing in x direction
dy = 0.01    # grid spacing in y direction
dz = 0.01    # grid spacing in z direction
dt = 0.001   # time step
T = 100      # simulation time
q = 1 # particle charge (in units of electron charge)
e = -q 
mu0 = 1
ep0 = 1
m = 1     # particle mass (in units of proton mass)
B0 = 10     # magnetic field strength
R = 4.0      # toroidal radius
a = 0.09      # plasma minor radius

q_over_m = q/m 
t = 0
tdata = []
xdata = []
ydata = []
zdata = []
plotRealTime = True


# Initialize particle positions, velocities, charges, and masses
x = random.uniform(0.01, a, Np)
y = random.uniform(0.01, a, Np)
z = random.uniform(0.01, a, Np)
vx = 1 * ones(Np)
vy = 1 * ones(Np)
vz = 1 * ones(Np)
Q = q * ones(Np)
m = m * ones(Np)

# Initialize electric and magnetic fields
Ex = zeros((int(Lx/dx), int(Ly/dy), int(Lz/dz)))
Ey = zeros((int(Lx/dx), int(Ly/dy), int(Lz/dz)))
Ez = zeros((int(Lx/dx), int(Ly/dy), int(Lz/dz)))
Bx = zeros((int(Lx/dx), int(Ly/dy), int(Lz/dz)))
By = zeros((int(Lx/dx), int(Ly/dy), int(Lz/dz)))
Bz = B0 * ones((int(Lx/dx), int(Ly/dy), int(Lz/dz)))






# Define functions for calculating charge density
def calculate_charge_density(x, y, z):
    rho = zeros((int(Lx/dx), int(Ly/dy), int(Lz/dz)))
    for i in range(Np):
    
        cn = [int(x[i]/dx),int(y[i]/dy),int(z[i]/dz)]
        
        n1 = (cn[0]-1,cn[1]-1,cn[2]-1)
        n2 = (cn[0]+1-1,cn[1]-1,cn[2]-1)
        n3 = (cn[0]+1-1,cn[1]+1-1,cn[2]-1)
        n4 = (cn[0]-1,cn[1]+1-1,cn[2]-1)
        n5 = (cn[0]-1,cn[1]-1,cn[2]-1+1)
        n6 = (cn[0]+1-1,cn[1]-1,cn[2]+1-1)
        n7 = (cn[0]+1-1,cn[1]+1-1,cn[2]+1-1)
        n8 = (cn[0]-1,cn[1]+1-1,cn[2]+1-1)
        
        w1 = abs((x[i]/dx-cn[0])*(y[i]/dy-cn[1])*(z[i]/dz-cn[2]))
        w2 = abs((x[i]/dx-cn[0]-1)*(y[i]/dy-cn[1])*(z[i]/dz-cn[2]))
        w3 = abs((x[i]/dx-cn[0]-1)*(y[i]/dy-cn[1]-1)*(z[i]/dz-cn[2]))
        w4 = abs((x[i]/dx-cn[0])*(y[i]/dy-cn[1]-1)*(z[i]/dz-cn[2]))
        w5 = abs((x[i]/dx-cn[0])*(y[i]/dy-cn[1])*(z[i]/dz-cn[2]-1))
        w6 = abs((x[i]/dx-cn[0]-1)*(y[i]/dy-cn[1])*(z[i]/dz-cn[2]-1))
        w7 = abs((x[i]/dx-cn[0]-1)*(y[i]/dy-cn[1]-1)*(z[i]/dz-cn[2]-1))
        w8 = abs((x[i]/dx-cn[0])*(y[i]/dy-cn[1]-1)*(z[i]/dz-cn[2]-1))
        
        rho[n1]+=Q[i]*(w7)
        rho[n2]+=Q[i]*(w8)
        rho[n3]+=Q[i]*(w5)
        rho[n4]+=Q[i]*(w6)
        rho[n5]+=Q[i]*(w3)
        rho[n6]+=Q[i]*(w4)
        rho[n7]+=Q[i]*(w1)
        rho[n8]+=Q[i]*(w2)
    return rho


def calculate_magnetic_field(Ex, Ey, Ez):
    for i in range(int(Lx/dx)):
        for j in range(int(Ly/dy)):
            for k in range (int(Lz/dz)):
                if (j+1 == Ly/dy) or (k+1 == Lz/dz):
                    Bx[i,j,k] = 0
                else:
                    Bx[i,j,k] += Ey[i,j,k+1] - Ey[i,j,k] - Ez[i,j+1,k] + Ez[i,j,k]
                if (i+1 == Lx/dx) or (k+1 == Lz/dz):      
                    By[i,j,k] = 0
                else:
                    By[i,j,k] += Ez[i+1,j,k] - Ez[i,j,k] - Ex[i,j,k+1] + Ex[i,j,k]
                if (i+1 == Lx/dx) or (j+1 == Ly/dy):
                    Bz[i,j,k] = 0
                else:
                    Bz[i,j,k] += Ex[i,j+1,k] - Ex[i,j,k] - Ey[i+1,j,k] + Ey[i,j,k]
    return Bx, By, Bz

#Calculating electric field using magnetic field and current density

def calculate_electric__field(Bx, By, Bz, Jx, Jy, Jz):
    for i in range(int(Lx/dx)):
        for j in range(int(Ly/dy)):
            for k in range (int(Lz/dz)):
                if (j-1 < 0) or (k-1 < 0):
                    Ex[i,j,k] = 0
                else:
                    Ex[i,j,k] += ((By[i,j,k-1] - By[i,j,k] - Bz[i,j-1,k] + Bz[i,j,k])/mu0 - Jx[i,j,k])/ep0
                if (i-1 < 0) or (k-1 < 0):      
                    Ey[i,j,k] = 0
                else:
                    Ey[i,j,k] += ((Bz[i-1,j,k] - Bz[i,j,k] - Bx[i,j,k-1] + Bx[i,j,k])/mu0 - Jy[i,j,k])/ep0
                if (i-1 < 0) or (j-1 < 0):
                    Ez[i,j,k] = 0
                else:
                    Ez[i,j,k] += ((Bx[i,j-1,k] - Bx[i,j,k] - By[i-1,j,k] + Ey[i,j,k])/mu0 - Jz[i,j,k])/ep0
    return Ex, Ey, Ez
    
#Calculating Current Density in the node using rho
def calculate_current_density(rho, x, y, z, vx, vy, vz):
    Jx = zeros((int(Lx/dx), int(Ly/dy), int(Lz/dz)))
    Jy = zeros((int(Lx/dx), int(Ly/dy), int(Lz/dz)))
    Jz = zeros((int(Lx/dx), int(Ly/dy), int(Lz/dz)))
    for i in range(Np):
    
        cn = [int(x[i]/dx),int(y[i]/dy),int(z[i]/dz)]
        
        n1 = (cn[0]-1,cn[1]-1,cn[2]-1)
        n2 = (cn[0]+1-1,cn[1]-1,cn[2]-1)
        n3 = (cn[0]+1-1,cn[1]+1-1,cn[2]-1)
        n4 = (cn[0]-1,cn[1]+1-1,cn[2]-1)
        n5 = (cn[0]-1,cn[1]-1,cn[2]-1+1)
        n6 = (cn[0]+1-1,cn[1]-1,cn[2]+1-1)
        n7 = (cn[0]+1-1,cn[1]+1-1,cn[2]+1-1)
        n8 = (cn[0]-1,cn[1]+1-1,cn[2]+1-1)
            
        Jx[n1] = rho[n1]*vx[i]
        Jx[n2] = rho[n2]*vx[i]
        Jx[n3] = rho[n3]*vx[i]
        Jx[n4] = rho[n4]*vx[i]
        Jx[n5] = rho[n5]*vx[i]
        Jx[n6] = rho[n6]*vx[i]
        Jx[n7] = rho[n7]*vx[i]
        Jx[n8] = rho[n8]*vx[i]
        
        Jy[n1] = rho[n1]*vy[i]
        Jy[n2] = rho[n2]*vy[i]
        Jy[n3] = rho[n3]*vy[i]
        Jy[n4] = rho[n4]*vy[i]
        Jy[n5] = rho[n5]*vy[i]
        Jy[n6] = rho[n6]*vy[i]
        Jy[n7] = rho[n7]*vy[i]
        Jy[n8] = rho[n8]*vy[i]
        
        Jz[n1] = rho[n1]*vz[i]
        Jz[n2] = rho[n2]*vz[i]
        Jz[n3] = rho[n3]*vz[i]
        Jz[n4] = rho[n4]*vz[i]
        Jz[n5] = rho[n5]*vz[i]
        Jz[n6] = rho[n6]*vz[i]
        Jz[n7] = rho[n7]*vz[i]
        Jz[n8] = rho[n8]*vz[i]
        
    return Jx, Jy, Jz      

# we now reverse interpolate the electric field on the particle position.
def reverse_interpolate_electric_field(Ex, Ey, Ez, x, y, z):
    #Displaced mesh form of Ex, Ey, Ez lattices    
    
    Ex_at_pos = []
        
    Ey_at_pos = []
    
    Ez_at_pos = []
        
    for i in range(Np):
    
        cn = [int(x[i]/dx),int(y[i]/dy),int(z[i]/dz)]
        #node points of Ex lattice wrt base lattice
        if (x[i]/dx - cn[0])>0.5:
            l1 = [cn[0]+0.5,cn[1],cn[2]]
            l2 = [cn[0]+1+0.5,cn[1],cn[2]]
            l3 = [cn[0]+1+0.5,cn[1]+1,cn[2]]
            l4 = [cn[0]+0.5,cn[1]+1,cn[2]]
            l5 = [cn[0]+0.5,cn[1],cn[2]+1]
            l6 = [cn[0]+1+0.5,cn[1],cn[2]+1]
            l7 = [cn[0]+1+0.5,cn[1]+1,cn[2]+1]
            l8 = [cn[0]+0.5,cn[1]+1,cn[2]+1]
        else:
            l1 = [cn[0]-0.5,cn[1],cn[2]]
            l2 = [cn[0]+1-0.5,cn[1],cn[2]]
            l3 = [cn[0]+1-0.5,cn[1]+1,cn[2]]
            l4 = [cn[0]-0.5,cn[1]+1,cn[2]]
            l5 = [cn[0]-0.5,cn[1],cn[2]+1]
            l6 = [cn[0]+1-0.5,cn[1],cn[2]+1]
            l7 = [cn[0]+1-0.5,cn[1]+1,cn[2]+1]
            l8 = [cn[0]-0.5,cn[1]+1,cn[2]+1]
            
            
        #node points of Ey lattice
        if (y[i]/dy - cn[1])>0.5:
            m1 = [cn[0],cn[1]+0.5,cn[2]]
            m2 = [cn[0]+1,cn[1]+0.5,cn[2]]
            m3 = [cn[0]+1,cn[1]+1+0.5,cn[2]]
            m4 = [cn[0],cn[1]+1+0.5,cn[2]]
            m5 = [cn[0],cn[1]+0.5,cn[2]+1]
            m6 = [cn[0]+1,cn[1]+0.5,cn[2]+1]
            m7 = [cn[0]+1,cn[1]+1+0.5,cn[2]+1]
            m8 = [cn[0],cn[1]+1+0.5,cn[2]+1]
        else:
            m1 = [cn[0],cn[1]-0.5,cn[2]]
            m2 = [cn[0]+1,cn[1]-0.5,cn[2]]
            m3 = [cn[0]+1,cn[1]+1-0.5,cn[2]]
            m4 = [cn[0],cn[1]+1-0.5,cn[2]]
            m5 = [cn[0],cn[1]-0.5,cn[2]+1]
            m6 = [cn[0]+1,cn[1]-0.5,cn[2]+1]
            m7 = [cn[0]+1,cn[1]+1-0.5,cn[2]+1]
            m8 = [cn[0],cn[1]+1-0.5,cn[2]+1]
            
            
        #node points of Ez lattice
        if (z[i]/dz - cn[2])>0.5:
            n1 = [cn[0],cn[1],cn[2]+0.5]
            n2 = [cn[0]+1,cn[1],cn[2]+0.5]
            n3 = [cn[0]+1,cn[1]+1,cn[2]+0.5]
            n4 = [cn[0],cn[1]+1,cn[2]+0.5]
            n5 = [cn[0],cn[1],cn[2]+1+0.5]
            n6 = [cn[0]+1,cn[1],cn[2]+1+0.5]
            n7 = [cn[0]+1,cn[1]+1,cn[2]+1+0.5]
            n8 = [cn[0],cn[1]+1,cn[2]+1+0.5]
        else:
            n1 = [cn[0],cn[1],cn[2]-0.5]
            n2 = [cn[0]+1,cn[1],cn[2]-0.5]
            n3 = [cn[0]+1,cn[1]+1,cn[2]-0.5]
            n4 = [cn[0],cn[1]+1,cn[2]-0.5]
            n5 = [cn[0],cn[1],cn[2]+1-0.5]
            n6 = [cn[0]+1,cn[1],cn[2]+1-0.5]
            n7 = [cn[0]+1,cn[1]+1,cn[2]+1-0.5]
            n8 = [cn[0],cn[1]+1,cn[2]+1-0.5]
            
            
        # this is for updating the cell number if the particle is found in a lower cell   
        if (x[i]/dx - cn[0])<0.5:
            cn[0] += -1
        if (y[i]/dy - cn[1])<0.5:
            cn[1] += -1
        if (z[i]/dz - cn[2])<0.5:
            cn[2] += -1
            
            
        pos = [x[i]/dx,y[i]/dy,z[i]/dz]
        
        #weights for Ex lattice
        wl1 = abs(prod(subtract(pos , l1)))
        wl2 = abs(prod(subtract(pos , l2)))
        wl3 = abs(prod(subtract(pos , l3)))
        wl4 = abs(prod(subtract(pos , l4)))
        wl5 = abs(prod(subtract(pos , l5)))
        wl6 = abs(prod(subtract(pos , l6)))
        wl7 = abs(prod(subtract(pos , l7)))
        wl8 = abs(prod(subtract(pos , l8)))
        # print(wl1+wl2+wl3+wl4+wl5+wl6+wl7+wl8)
        
        #weights for Ey lattice
        wm1 = abs(prod(subtract(pos , m1)))
        wm2 = abs(prod(subtract(pos , m2)))
        wm3 = abs(prod(subtract(pos , m3)))
        wm4 = abs(prod(subtract(pos , m4)))
        wm5 = abs(prod(subtract(pos , m5)))
        wm6 = abs(prod(subtract(pos , m6)))
        wm7 = abs(prod(subtract(pos , m7)))
        wm8 = abs(prod(subtract(pos , m8)))
        # print(wm1+wm2+wm3+wm4+wm5+wm6+wm7+wm8)
        
        #weights for Ez lattice
        wn1 = abs(prod(subtract(pos , n1)))
        wn2 = abs(prod(subtract(pos , n2)))
        wn3 = abs(prod(subtract(pos , n3)))
        wn4 = abs(prod(subtract(pos , n4)))
        wn5 = abs(prod(subtract(pos , n5)))
        wn6 = abs(prod(subtract(pos , n6)))
        wn7 = abs(prod(subtract(pos , n7)))
        wn8 = abs(prod(subtract(pos , n8)))
        # print(wn1+wn2+wn3+wn4+wn5+wn6+wn7+wn8)
        
    #actual node points of Ex, Ey, Ez lattices
        a1 = (cn[0]-1,cn[1]-1,cn[2]-1)
        a2 = (cn[0]+1-1,cn[1]-1,cn[2]-1)
        a3 = (cn[0]+1-1,cn[1]+1-1,cn[2]-1)
        a4 = (cn[0]-1,cn[1]+1-1,cn[2]-1)
        a5 = (cn[0]-1,cn[1]-1,cn[2]+1-1)
        a6 = (cn[0]+1-1,cn[1]-1,cn[2]+1-1)
        a7 = (cn[0]+1-1,cn[1]+1-1,cn[2]+1-1)
        a8 = (cn[0]-1,cn[1]+1-1,cn[2]+1-1)
        
        #electric field at the position of the particle
        
        Ex_at_pos.append(Ex[a1]*wl7 + Ex[a2]*wl8 + Ex[a3]*wl5 + Ex[a4]*wl6 + Ex[a5]*wl3 + Ex[a6]*wl4 + Ex[a7]*wl1 + Ex[a8]*wl2)
        
        Ey_at_pos.append(Ey[a1]*wm7 + Ey[a2]*wm8 + Ey[a3]*wm5 + Ey[a4]*wm6 + Ey[a5]*wm3 + Ey[a6]*wm4 + Ey[a7]*wm1 + Ey[a8]*wm2)
        
        Ez_at_pos.append(Ez[a1]*wn7 + Ez[a2]*wn8 + Ez[a3]*wn5 + Ez[a4]*wn6 + Ez[a5]*wn3 + Ez[a6]*wn4 + Ez[a7]*wn1 + Ez[a8]*wn2)
        
        #print(Ex_at_pos, Ey_at_pos, Ez_at_pos)
        
    return Ex_at_pos, Ey_at_pos, Ez_at_pos
        

# we now reverse interpolate the magnetic field on the particle position.
def reverse_interpolate_magnetic_field(Bx, By, Bz, x, y, z):
    #Displaced mesh form of Bx, By, Bz lattices    
    
    Bx_at_pos = []
    
    By_at_pos = []
    
    Bz_at_pos = []
        
    for i in range(Np):
    
        cn = [int(x[i]/dx),int(y[i]/dy),int(z[i]/dz)]
        #node points of Bx lattice wrt base lattice
        if (z[i]/dz - cn[2])>0.5 and (y[i]/dy - cn[1])>0.5:
            l1 = [cn[0],cn[1]+0.5,cn[2]+0.5]
            l2 = [cn[0]+1,cn[1]+0.5,cn[2]+0.5]
            l3 = [cn[0]+1,cn[1]+1+0.5,cn[2]+0.5]
            l4 = [cn[0],cn[1]+1+0.5,cn[2]+0.5]
            l5 = [cn[0],cn[1]+0.5,cn[2]+1+0.5]
            l6 = [cn[0]+1,cn[1]+0.5,cn[2]+1+0.5]
            l7 = [cn[0]+1,cn[1]+1+0.5,cn[2]+1+0.5]
            l8 = [cn[0],cn[1]+1+0.5,cn[2]+1+0.5]
            
        elif (z[i]/dz - cn[2])>0.5 and (y[i]/dy - cn[1])<0.5:
            l1 = [cn[0],cn[1]-0.5,cn[2]+0.5]
            l2 = [cn[0]+1,cn[1]-0.5,cn[2]+0.5]
            l3 = [cn[0]+1,cn[1]+1-0.5,cn[2]+0.5]
            l4 = [cn[0],cn[1]+1-0.5,cn[2]+0.5]
            l5 = [cn[0],cn[1]-0.5,cn[2]+1+0.5]
            l6 = [cn[0]+1,cn[1]-0.5,cn[2]+1+0.5]
            l7 = [cn[0]+1,cn[1]+1-0.5,cn[2]+1+0.5]
            l8 = [cn[0],cn[1]+1-0.5,cn[2]+1+0.5]
            
        elif (z[i]/dz - cn[2])<0.5 and (y[i]/dy - cn[1])>0.5:
            l1 = [cn[0],cn[1]+0.5,cn[2]-0.5]
            l2 = [cn[0]+1,cn[1]+0.5,cn[2]-0.5]
            l3 = [cn[0]+1,cn[1]+1+0.5,cn[2]-0.5]
            l4 = [cn[0],cn[1]+1+0.5,cn[2]-0.5]
            l5 = [cn[0],cn[1]+0.5,cn[2]+1-0.5]
            l6 = [cn[0]+1,cn[1]+0.5,cn[2]+1-0.5]
            l7 = [cn[0]+1,cn[1]+1+0.5,cn[2]+1-0.5]
            l8 = [cn[0],cn[1]+1+0.5,cn[2]+1-0.5]
            
        elif (z[i]/dz - cn[2])<0.5 and (y[i]/dy - cn[1])<0.5:
            l1 = [cn[0],cn[1]-0.5,cn[2]-0.5]
            l2 = [cn[0]+1,cn[1]-0.5,cn[2]-0.5]
            l3 = [cn[0]+1,cn[1]+1-0.5,cn[2]-0.5]
            l4 = [cn[0],cn[1]+1-0.5,cn[2]-0.5]
            l5 = [cn[0],cn[1]-0.5,cn[2]+1-0.5]
            l6 = [cn[0]+1,cn[1]-0.5,cn[2]+1-0.5]
            l7 = [cn[0]+1,cn[1]+1-0.5,cn[2]+1-0.5]
            l8 = [cn[0],cn[1]+1-0.5,cn[2]+1-0.5]
            
            
        #node points of By lattice wrt base lattice
        if (z[i]/dz - cn[2])>0.5 and (x[i]/dx - cn[0])>0.5:
            m1 = [cn[0]+0.5,cn[1],cn[2]+0.5]
            m2 = [cn[0]+1+0.5,cn[1],cn[2]+0.5]
            m3 = [cn[0]+1+0.5,cn[1]+1,cn[2]+0.5]
            m4 = [cn[0]+0.5,cn[1]+1,cn[2]+0.5]
            m5 = [cn[0]+0.5,cn[1],cn[2]+1+0.5]
            m6 = [cn[0]+1+0.5,cn[1],cn[2]+1+0.5]
            m7 = [cn[0]+1+0.5,cn[1]+1,cn[2]+1+0.5]
            m8 = [cn[0]+0.5,cn[1]+1,cn[2]+1+0.5]
            
        elif (z[i]/dz - cn[2])>0.5 and (x[i]/dx - cn[0])<0.5:
            m1 = [cn[0]-0.5,cn[1],cn[2]+0.5]
            m2 = [cn[0]+1-0.5,cn[1],cn[2]+0.5]
            m3 = [cn[0]+1-0.5,cn[1]+1,cn[2]+0.5]
            m4 = [cn[0]-0.5,cn[1]+1,cn[2]+0.5]
            m5 = [cn[0]-0.5,cn[1],cn[2]+1+0.5]
            m6 = [cn[0]+1-0.5,cn[1],cn[2]+1+0.5]
            m7 = [cn[0]+1-0.5,cn[1]+1,cn[2]+1+0.5]
            m8 = [cn[0]-0.5,cn[1]+1,cn[2]+1+0.5]
            
        elif (z[i]/dz - cn[2])<0.5 and (x[i]/dx - cn[0])>0.5:
            m1 = [cn[0]+0.5,cn[1],cn[2]-0.5]
            m2 = [cn[0]+1+0.5,cn[1],cn[2]-0.5]
            m3 = [cn[0]+1+0.5,cn[1]+1,cn[2]-0.5]
            m4 = [cn[0]+0.5,cn[1]+1,cn[2]-0.5]
            m5 = [cn[0]+0.5,cn[1],cn[2]+1-0.5]
            m6 = [cn[0]+1+0.5,cn[1],cn[2]+1-0.5]
            m7 = [cn[0]+1+0.5,cn[1]+1,cn[2]+1-0.5]
            m8 = [cn[0]+0.5,cn[1]+1,cn[2]+1-0.5]
            
        elif (z[i]/dz - cn[2])<0.5 and (x[i]/dx - cn[0])<0.5:
            m1 = [cn[0]-0.5,cn[1],cn[2]-0.5]
            m2 = [cn[0]+1-0.5,cn[1],cn[2]-0.5]
            m3 = [cn[0]+1-0.5,cn[1]+1,cn[2]-0.5]
            m4 = [cn[0]-0.5,cn[1]+1,cn[2]-0.5]
            m5 = [cn[0]-0.5,cn[1],cn[2]+1-0.5]
            m6 = [cn[0]+1-0.5,cn[1],cn[2]+1-0.5]
            m7 = [cn[0]+1-0.5,cn[1]+1,cn[2]+1-0.5]
            m8 = [cn[0]-0.5,cn[1]+1,cn[2]+1-0.5]
            
        
            
        #node points of Bz lattice wrt base lattice
        if (y[i]/dy - cn[1])>0.5 and (x[i]/dx - cn[0])>0.5:
            n1 = [cn[0]+0.5,cn[1]+0.5,cn[2]]
            n2 = [cn[0]+1+0.5,cn[1]+0.5,cn[2]]
            n3 = [cn[0]+1+0.5,cn[1]+1+0.5,cn[2]]
            n4 = [cn[0]+0.5,cn[1]+1+0.5,cn[2]]
            n5 = [cn[0]+0.5,cn[1]+0.5,cn[2]+1]
            n6 = [cn[0]+1+0.5,cn[1]+0.5,cn[2]+1]
            n7 = [cn[0]+1+0.5,cn[1]+1+0.5,cn[2]+1]
            n8 = [cn[0]+0.5,cn[1]+1+0.5,cn[2]+1]
            
        elif (y[i]/dy - cn[1])>0.5 and (x[i]/dx - cn[0])<0.5:
            n1 = [cn[0]-0.5,cn[1]+0.5,cn[2]]
            n2 = [cn[0]+1-0.5,cn[1]+0.5,cn[2]]
            n3 = [cn[0]+1-0.5,cn[1]+1+0.5,cn[2]]
            n4 = [cn[0]-0.5,cn[1]+1+0.5,cn[2]]
            n5 = [cn[0]-0.5,cn[1]+0.5,cn[2]+1]
            n6 = [cn[0]+1-0.5,cn[1]+0.5,cn[2]+1]
            n7 = [cn[0]+1-0.5,cn[1]+1+0.5,cn[2]+1]
            n8 = [cn[0]-0.5,cn[1]+1+0.5,cn[2]+1]
            
        elif (y[i]/dy - cn[1])<0.5 and (x[i]/dx - cn[0])>0.5:
            n1 = [cn[0]+0.5,cn[1]-0.5,cn[2]]
            n2 = [cn[0]+1+0.5,cn[1]-0.5,cn[2]]
            n3 = [cn[0]+1+0.5,cn[1]+1-0.5,cn[2]]
            n4 = [cn[0]+0.5,cn[1]+1-0.5,cn[2]]
            n5 = [cn[0]+0.5,cn[1]-0.5,cn[2]+1]
            n6 = [cn[0]+1+0.5,cn[1]-0.5,cn[2]+1]
            n7 = [cn[0]+1+0.5,cn[1]+1-0.5,cn[2]+1]
            n8 = [cn[0]+0.5,cn[1]+1-0.5,cn[2]+1]
            
        elif (y[i]/dy - cn[1])<0.5 and (x[i]/dx - cn[0])<0.5:
            n1 = [cn[0]-0.5,cn[1]-0.5,cn[2]]
            n2 = [cn[0]+1-0.5,cn[1]-0.5,cn[2]]
            n3 = [cn[0]+1-0.5,cn[1]+1-0.5,cn[2]]
            n4 = [cn[0]-0.5,cn[1]+1-0.5,cn[2]]
            n5 = [cn[0]-0.5,cn[1]-0.5,cn[2]+1]
            n6 = [cn[0]+1-0.5,cn[1]-0.5,cn[2]+1]
            n7 = [cn[0]+1-0.5,cn[1]+1-0.5,cn[2]+1]
            n8 = [cn[0]-0.5,cn[1]+1-0.5,cn[2]+1]
            
        # this is for updating the cell number if the particle is found in a lower cell    
        if (x[i]/dx - cn[0])<0.5:
            cn[0] += -1
        if (y[i]/dy - cn[1])<0.5:
            cn[1] += -1
        if (z[i]/dz - cn[2])<0.5:
            cn[2] += -1
            
            
        pos = [x[i]/dx,y[i]/dy,z[i]/dz]
        
        #weights for Bx lattice
        wl1 = abs(prod(subtract(pos , l1)))
        wl2 = abs(prod(subtract(pos , l2)))
        wl3 = abs(prod(subtract(pos , l3)))
        wl4 = abs(prod(subtract(pos , l4)))
        wl5 = abs(prod(subtract(pos , l5)))
        wl6 = abs(prod(subtract(pos , l6)))
        wl7 = abs(prod(subtract(pos , l7)))
        wl8 = abs(prod(subtract(pos , l8)))
        # print(wl1+wl2+wl3+wl4+wl5+wl6+wl7+wl8)
        
        #weights for By lattice
        wm1 = abs(prod(subtract(pos , m1)))
        wm2 = abs(prod(subtract(pos , m2)))
        wm3 = abs(prod(subtract(pos , m3)))
        wm4 = abs(prod(subtract(pos , m4)))
        wm5 = abs(prod(subtract(pos , m5)))
        wm6 = abs(prod(subtract(pos , m6)))
        wm7 = abs(prod(subtract(pos , m7)))
        wm8 = abs(prod(subtract(pos , m8)))
        # print(wm1+wm2+wm3+wm4+wm5+wm6+wm7+wm8)
        
        #weights for Bz lattice
        wn1 = abs(prod(subtract(pos , n1)))
        wn2 = abs(prod(subtract(pos , n2)))
        wn3 = abs(prod(subtract(pos , n3)))
        wn4 = abs(prod(subtract(pos , n4)))
        wn5 = abs(prod(subtract(pos , n5)))
        wn6 = abs(prod(subtract(pos , n6)))
        wn7 = abs(prod(subtract(pos , n7)))
        wn8 = abs(prod(subtract(pos , n8)))
        # print(wn1+wn2+wn3+wn4+wn5+wn6+wn7+wn8)
        
    #actual node points of Bx, By, Bz lattices
        a1 = (cn[0]-1,cn[1]-1,cn[2]-1)
        a2 = (cn[0]+1-1,cn[1]-1,cn[2]-1)
        a3 = (cn[0]+1-1,cn[1]+1-1,cn[2]-1)
        a4 = (cn[0]-1,cn[1]+1-1,cn[2]-1)
        a5 = (cn[0]-1,cn[1]-1,cn[2]+1-1)
        a6 = (cn[0]+1-1,cn[1]-1,cn[2]+1-1)
        a7 = (cn[0]+1-1,cn[1]+1-1,cn[2]+1-1)
        a8 = (cn[0]-1,cn[1]+1-1,cn[2]+1-1)
        
        #Magnetic field at the position of the particle
        
        Bx_at_pos.append(Bx[a1]*wl7 + Bx[a2]*wl8 + Bx[a3]*wl5 + Bx[a4]*wl6 + Bx[a5]*wl3 + Bx[a6]*wl4 + Bx[a7]*wl1 + Bx[a8]*wl2)
        
        By_at_pos.append(By[a1]*wm7 + By[a2]*wm8 + By[a3]*wm5 + By[a4]*wm6 + By[a5]*wm3 + By[a6]*wm4 + By[a7]*wm1 + By[a8]*wm2)
        
        Bz_at_pos.append(Bz[a1]*wn7 + Bz[a2]*wn8 + Bz[a3]*wn5 + Bz[a4]*wn6 + Bz[a5]*wn3 + Bz[a6]*wn4 + Bz[a7]*wn1 + Bz[a8]*wn2)
        
        # print(Bx_at_pos, By_at_pos, Bz_at_pos)
        
    return Bx_at_pos, By_at_pos, Bz_at_pos   

def position_vector_creator(x, y, z):  
    posi = [x[0], y[0], z[0]] 
    for i in range(1, len(x)):
        pp = [x[i], y[i], z[i]]
        posi = vstack((posi, pp))
        
    return posi

def velocity_vector_creator(vx, vy, vz):  
    v = [vx[0], vy[0], vz[0]] 
    for i in range(1, len(vx)):
        vv = [vx[i], vy[i], vz[i]]
        v = vstack((v, vv))
        
    return v

def electric_field_vector_creator(Ex_at_pos, Ey_at_pos, Ez_at_pos):  
    E = [Ex_at_pos[0], Ey_at_pos[0], Ez_at_pos[0]]   
    for i in range(1, len(Ex_at_pos)):
        e = [Ex_at_pos[i], Ey_at_pos[i], Ez_at_pos[i]]
        E = vstack((E,e))
        
    return E

def magnetic_field_vector_creator(Bx_at_pos, By_at_pos, Bz_at_pos):
    B = [Bx_at_pos[0],By_at_pos[0],Bz_at_pos[0]]
    for i in range(1, len(Bx_at_pos)):
        b = [Bx_at_pos[i],By_at_pos[i],Bz_at_pos[i]]
        B = vstack((B,b))     
        
    return B   
        
        

# Main Loop
for i in range(10):
    
    # this calculate the charge density, current density, magnetic field and electric field in the mesh
    rho = calculate_charge_density(x,y,z)
    Jx, Jy, Jz = calculate_current_density(rho, x, y, z, vx, vy, vz)
    Bx, By, Bz = calculate_magnetic_field(Ex, Ey, Ez)
    Ex, Ey, Ez = calculate_electric__field(Bx, By, Bz, Jx, Jy, Jz)
    
    #this interpolates the electric and magnetic fields at the position of the particles
    
    Ex_at_pos, Ey_at_pos, Ez_at_pos = reverse_interpolate_electric_field(Ex, Ey, Ez, x, y, z)
    Bx_at_pos, By_at_pos, Bz_at_pos = reverse_interpolate_magnetic_field(Bx, By, Bz, x, y, z)
    
    E = electric_field_vector_creator(Ex_at_pos, Ey_at_pos, Ez_at_pos)
    B = magnetic_field_vector_creator(Bx_at_pos, By_at_pos, Bz_at_pos)
    
    v = velocity_vector_creator(vx, vy, vz)
    v += q_over_m * E*dt/4   

    v_prime = v + q_over_m*cross(v, B)*dt/2
    vnew = v_prime + q_over_m * E * dt/4
    pos = position_vector_creator(x,y,z)
    pos += vnew*dt
    
    for i in range(Np):
        for j in range(3):
            if pos[i,j] > Lx:
                if pos[i,j] - Lx < 0.01:
                    pos[i,j] = 2*Lx -pos[i,j]
                    vnew[i,j] = -vnew[i,j]
                else:
                    pos[i,j] = 0.05
            if pos[i,j] < 0:
                if abs(pos[i,j]) < 0.01:
                    pos[i,j] = abs(pos[i,j])
                    vnew[i,j] = -vnew[i,j]
                else:
                    pos[i,j] = 0.05
    
    x = pos[:,0]
    y = pos[:,1]
    z = pos[:,2]
    
    vx = vnew[:,0]
    vy = vnew[:,1]
    vz = vnew[:,2]
    
    #record data at each step
    tdata.append(t)
    xdata.append((pos[:,0]))
    ydata.append((pos[:,1]))
    zdata.append((pos[:,2]))
    


#for animation
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')


def update(frame):
    ax.clear()
    ax.scatter(xdata[:frame], ydata[:frame], zdata[:frame])
    ax.set_xlim([0, Lx])
    ax.set_ylim([0, Ly])
    ax.set_zlim([0, Lz])
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('Particle Trajectory')


ani = FuncAnimation(fig, update, frames=len(x), interval=50)

ani.save('animation1.mp4', writer='ffmpeg')

plt.show()


    

    
    
    
