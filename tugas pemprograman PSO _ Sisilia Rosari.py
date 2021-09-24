import random
import numpy as np
import pandas as pd
from collections import Counter

# parameters
n_particles = 20
n_iteration = 200
numofcar = 3
capofcar= 100
numOfDest = 15

W  = 0.7
c1 = 1
c2 = 2

destination = np.arange(1,numOfDest+1)

particles = np.zeros((n_particles, numOfDest), dtype=int) 
velocities = np.zeros((n_particles, 15))

for i in range(n_particles):
    particles[i,:] = np.random.choice(destination, numOfDest, replace=False)
for i in range(n_particles):
    velocities[i,:] = np.random.uniform(-1,1,15)

#demand each route in particle
def getdemand(particles):
    demandpar = np.zeros((numOfDest), dtype=int)
    for j in range(numOfDest):
        demandpar[j]= demand[particles[j]-1][0]
    return demandpar

#constraint three car with capacity constraint 100 kg
def car(particles, demandpar):
    car_particle = np.zeros(numOfDest, dtype=int)
    demandcopy = demandpar.copy()
    for k in range(numofcar):
            mobil = capofcar        
            for j in range(numOfDest):  
                if demandcopy[j] < mobil:
                    mobil = mobil - demandcopy[j]
                    car_particle[j] = k
                    demandcopy[j] = 9999
                else:
                   continue
    demandcopy = demandpar.copy()            
    return car_particle

def diskrit(new_position):    
    interval_x = (max(new_position)-min(new_position))/new_position.shape[0]
    interval = []
    for i in range(new_position.shape[0]):
        a = min(new_position)+i*interval_x
        interval.append(a)
    interval = np.array(interval)        
    index=[] 
    for i in range(new_position.shape[0]):
        for j in range(new_position.shape[0]):
            if new_position[i] < interval[j]:
                break 
        index.append(j+1)
    index = np.array(index)
    return index

def objective_function(particles, car_particle):
    fitness=0
    for i in range(numofcar):
        no_mobil = np.where(car_particle==i)[0]
        totaljarak = dist[0,particles[no_mobil[0]]]
        for j in range(len(no_mobil)-1):
            start = particles[no_mobil[j]]
            end = particles[no_mobil[j]+1]
            totaljarak = totaljarak + dist[start,end]
        totaljarak = totaljarak + dist[end,0]
        fitness= fitness+totaljarak
    return int(fitness)

# check boundary violation
def boundary_violation(index):
    frequency = [item for item, count in Counter(index).items() if count > 1]
    missing = destination[~np.isin(destination, index)]
    index_frequency = []
    for i in range(len(frequency)):
        z = np.where(index==frequency[i])[0]
        for j in range(len(z)):
            if j>0:
                index_frequency.append(z[j])
    index_frequency.sort()
    for i in range(len(missing)):
        index[index_frequency[i]] = missing[i]
    return index

# Getting distances data 
dist = pd.read_csv('rutejarak.csv', header=0, index_col=0)
dist = dist.values

# Getting demand data
demand = pd.read_csv('demand.csv', header=0, index_col=0)
demand = demand.values

# for iteration 
pbest_location = particles.copy()
pbest_fitness_value = np.array([float('inf') for _ in range(n_particles)])
gbest_fitness_value = float('inf')
gbest_location = np.array([float('inf'), float('inf')])

gbest_history = np.zeros((n_iteration, 1))
gbest_location_history = np.zeros((n_iteration,15))

it = 0
# while maximum iteration or convergence criteria is not met
while it < n_iteration:	
    print('------------------------ ')
    print('ITERASI #' +str(it))
    print('>>>>> ')
    # For each particle:
    # Take fitness value Fitness(i,t) corresponded to location particle(i,t)
    for i in range (n_particles): 	
        demandpar = getdemand(particles[i])
        car_particle = car(particles[i],demandpar)
        fitness_value = objective_function(particles[i],car_particle)
        print(particles[i], ' ', fitness_value)
        
        if(pbest_fitness_value[i] > fitness_value):
            pbest_fitness_value[i] = fitness_value
            pbest_location[i] = particles[i].copy()
            
        if(gbest_fitness_value > fitness_value):
            gbest_fitness_value = fitness_value
            gbest_location = particles[i].copy()
    
    # Record gbest history 
    gbest_history[it,0] = gbest_fitness_value
    gbest_location_history[it,:] = gbest_location
    print(' ') 
    print('gbest')
    print(gbest_location, ' ', gbest_fitness_value)
    print(' ') 
    print('gbest demand')
    print(getdemand(gbest_location))
    print(' ') 
    print('gbest kendaraan')
    print(car(gbest_location, getdemand(gbest_location)))
            
    # For each particle:
    # Calculate particle velocity
    # Update all particle    
    for i in range(n_particles):
        r1 = round(random.random(),4)
        r2 = round(random.random(),4) 
        new_velocity = (W*velocities[i]) + (c1*r1*(pbest_location[i] - particles[i])) +                        (c2*r2*(gbest_location-particles[i])) 
        new_position = new_velocity.round(4) + particles[i]
        new_position = diskrit(new_position)
        
        #new_position = checkfeasibility(new_position)
        #check boundary violation
        new_position = boundary_violation(new_position)
        particles[i] = new_position.copy()
        velocities[i] = new_velocity.copy()
        demandpar = getdemand(particles)
        car_particle = car(particles, demandpar)
    it += 1
    print(' ') #end while

# gbest
print('---------------------------------------')
print('Final Gbest')
print(gbest_location, ' ', gbest_fitness_value)    
plt = pd.DataFrame(gbest_history, columns=['Gbest']).plot()
demandpar = getdemand(gbest_location)
carpar = car(gbest_location, demandpar)
for i in range(numofcar):
    no_mobil = np.where(carpar==i)[0]
    print('rute dari mobil', i+1, 'adalah', gbest_location[no_mobil])


