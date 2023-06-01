from numpy import *
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation



for po in range(3):
    start_time = time.time()
    np = 10**(po+1)
    delta_t = 0.1
    dimensions = [1000, 1000, 1000]
    charge = random.randint(-1,1)

    def vector_magnitude(vector):
        #Cacluclates the magnitude of the vector
        mag = 0.0
        for component in vector:
            mag = mag + (component ** 2)
        mag = mag ** 0.5
        return mag

    def unit_vector(vector):
        if vector_magnitude(vector) != 0:
            u_vector = vector/vector_magnitude(vector)
        else:
            u_vector = np.zeroes(vector.shape)
        return u_vector

    def electrostatic(pos1,pos2):
        
        distance =  pos1 - pos2
        force = 9e9 * charge / (vector_magnitude(distance)**2)
        force = force*unit_vector(distance)
        return force




    def update():
        force = electrostatic(pos1, pos2)
        acc = force/(mass)
        vel += acc*delta_t
        pos1 = pos1 + vel*delta_t
        collide()
        return pos1
        
    def collide(pos1, veli):
        if pos1[0]<= 0 or pos1[0]>=dimensions[0]:
            veli *= asarray([-1,1,1])
        if pos1[1] <=0 or pos1[1] >= dimensions[1]:
            veli *= asarray([1,-1,1])
        if pos1[2]<= 0 or pos1[2]>=dimensions[2]:
            veli *= asarray([1,1,-1])
        return veli 
        
    def velocity_vector_creator(vx, vy, vz):  
        v = [vx[0], vy[0], vz[0]] 
        for i in range(1, len(vx)):
            vv = [vx[i], vy[i], vz[i]]
            v = vstack((v, vv))
            
        return v
            
    x = random.uniform(0.01, dimensions[0], np)
    y = random.uniform(0.01, dimensions[1], np)
    z = random.uniform(0.01, dimensions[2], np)
    vx = random.uniform(-5, 5, np)
    vy = random.uniform(-5, 5, np)
    vz = random.uniform(-5, 5, np)
    mass = 1
    force = 0
    vel = velocity_vector_creator(vx, vy, vz)
    for t in range(10):
        for i in range(np):
            for j in range(np):
                if i == j:
                    pass
                else:  
                    pos1 = asarray([x[i],y[i],z[i]])
                    pos2 = asarray([x[j],y[j],z[j]])
                    # print(pos1)
                    force += electrostatic(pos1, pos2)
            acc = force/(mass)
            veli = vel[i,:]
            veli += acc*delta_t
            pos1 = pos1 + veli*delta_t
            #print(pos1)
            collide(pos1, veli)
            vel[i,:] = veli
            [x[i],y[i],z[i]] = pos1
    
    
            
         
            
            
        
