# -*- coding: utf-8 -*-
"""
MVP Checkpoint 1

@author: Surface
"""
#imports 
import sys
import numpy as np 
import math
import random 
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def pbc(pos, l):
    # Use mod function to find position of particle within pbc
    return np.mod(pos, l)

#calculate system energy due to nearest neighbours 
def energy(lattice, x, y,):
    centre = lattice[x,y]
    up = lattice[x,pbc(y+1, Ny)]
    down = lattice[x,pbc(y-1, Ny)]
    left = lattice[pbc(x-1, Nx),y]
    right = lattice[pbc(x+1, Nx),y]
    energy = -1*(centre*(up+down+left+right))
    return energy 

#constants
J = 1
k_b = 1

steps = int(input(">>>number of steps = ")) 

dynamic = input(">>>Glauber (G) or Kawasaki (K) = ") 

#define lattice dimenasion 
N = int(input(">>>lattice dimensions = ")) 
Nx = N
Ny = N

Ti = float(input(">>>initial tmeperature = "))
Tf = float(input(">>>final temperature = "))
t = float(input(">>>temperature interval = "))

animate = input(">>>animate (Y/N) = ")

#create empty lattice 
#lattice = np.zeros((Nx,Ny), dtype=int)

#iterate through points of the lattice 
#for i in range(Nx):
#    for j in range(Ny):
#        #generate a random number between 0 and 1
#        r = random.random()
#        #if r<0.5 assign lattice point value -1
#        if (r<0.5):
#            lattice[i,j] = -1
#        #if r>0.5 assign lattice point value 1
#        if (r>0.5):
#            lattice[i,j] = 1



#show lattice
#fig = plt.figure()
#im = plt.imshow(lattice, animated = True)
 
#initiate variable lists
avg_mag_lst = []
suscept_lst = []
avg_energy_lst = []
avg_energy_sq_lst = []
heat_cap_lst = [] 
temp = [] 
errors = []

#Glauber dynamics 
if dynamic == "G":
    #set to ground state configuration (spins all one way, e.g. up)
    lattice = np.ones((Nx,Ny), dtype=int)
    print(lattice)
    #initial temperature
    T = Ti
    #while loop for temperature range
    while T <= Tf:  
        #append temperature to list
        temp.append(T)
        print(T)
        #initiate magnetisation and energy lists
        mag = []
        energy_all = []
        heat_cap_s_lst = []
        #for number of steps
        for n in range(steps):
            #one sweep is NxN measurements 
            for a in range(Nx*Ny):
                
                #generate random values for coordinates of spin reverse
                x = random.randint(0, Nx-1)
                y = random.randint(0, Ny-1)
                
                #set new spin of random coordinate
                new_spin = -1*lattice[x,y]
                    
                #print([x,y])
                
                #calculate change in energy of the system 
                delta_energy = -2*energy(lattice, x, y)
                                
                #metrapolis algorithm - determine when to allow spin reverse
                if delta_energy < 0:
                    lattice[x,y] = new_spin
                else:
                    prob = math.exp(-delta_energy/T)
                    if random.random() <= prob:
                        lattice[x,y] = new_spin
                
            #calculate magnetisation and energy every 10 steps    
            if n >= 100 and n % 10 == 0:
                mag.append(np.sum(lattice))
                energy_tot = 0
                for i in range(Nx):
                    for j in range(Ny):
                        energy_tot += energy(lattice, i, j)
                energy_all.append(0.5*energy_tot)
           
            if animate == 'Y':        
                if(n % 100 == 0):        
                   f = open('spins_g.dat','w')
                   for i in range(Nx):
                       for j in range(Ny):
                           f.write('%d %d %lf\n'%(i,j,lattice[i,j]))
                   f.close()
                   #show animation
                   plt.cla()
                   im=plt.imshow(lattice, animated=True)
                   plt.draw()
                   plt.pause(0.0001) 
        

         
        #calculate average magnetisation and suceptibility across temperatures
        mag_sq = np.array(mag)**2
        avg_mag = np.mean(np.abs(mag))
        avg_mag_lst.append(avg_mag)
        suscept = (np.mean(mag_sq)-(np.mean(mag)**2))/(T*(N**2))
        suscept_lst.append(suscept)
        
        #calculate average energy and heat capacity across temperatures
        energy_all_sq = np.array(energy_all)**2
        avg_energy = np.mean(energy_all)
        avg_energy_lst.append(avg_energy)
        heat_cap = (np.mean(energy_all_sq)-(avg_energy**2))/((N**2)*(T**2))
        heat_cap_lst.append(heat_cap)
        
        #bootstrap error
        m = 0
        heat_cap_s_lst = []
        while m <= 1000:
            #get sample (with replacment)
            s = np.random.choice(energy_all, len(energy_all))
            avg_s = np.mean(s)
            avg_s_sq = np.array(s)**2
            #calculate heat capacity for sample
            heat_cap_s = (np.mean(avg_s_sq)-(avg_s**2))/((N**2)*(T**2))
            #append heat capacity to a list of all sample heat capacities
            heat_cap_s_lst.append(heat_cap_s)
            m += 1
            #calculate variance and standard deviation    
        var = (np.mean(np.array(heat_cap_s_lst)**2)-np.mean(heat_cap_s_lst)**2)
        std = var**(1/2) 
        errors.append(std)
                    
        #increase the temperature
        T += t
    
    print(errors)
        
    #calculate variance and standard deviation    
    var = (np.mean(np.array(heat_cap_s_lst)**2)-np.mean(heat_cap_s_lst)**2)
    std = var**(1/2)    
    
    #write to files
    f = open('energy_g.dat','w')
    for v in range(len(temp)):
            f.write(str(temp[v]) + ": " + str(avg_energy_lst[v]) + "\n")
    f.close()    
    
    f = open('magnetisation_g.dat','w')
    for v in range(len(temp)):
            f.write(str(temp[v]) + ": " + str(avg_mag_lst[v]) + "\n")
    f.close()   
    
    f = open('susceptibility_g.dat','w')
    for v in range(len(temp)):
            f.write(str(temp[v]) + ": " + str(suscept_lst[v]) + "\n")
    f.close()  
    
    f = open('heat_capacity_g.dat','w')
    for v in range(len(temp)):
            f.write(str(temp[v]) + ": " + str(heat_cap_lst[v]) + " +/- " + str(errors[v]) + "\n")
    f.close()
    
    #produce plots
    plt.title("Magnetisation change with temperature")
    plt.xlabel("Temperature")
    plt.ylabel("Magnetisation")
    plt.plot(temp, avg_mag_lst)
    plt.savefig('magnetisation_g.png')
    plt.show()

    
    plt.title("Energy change with temperature")
    plt.xlabel("Temperature")
    plt.ylabel("Energy")
    plt.plot(temp, avg_energy_lst)
    plt.savefig('energy_g.png')
    plt.show()

    
    plt.title("Susecptibility change with temperature")
    plt.xlabel("Temperature")
    plt.ylabel("Susceptibility")
    plt.plot(temp, suscept_lst)
    plt.savefig('susceptibility_g.png')
    plt.show()

    
    plt.title("Heat Capacity change with temperature")
    plt.xlabel("Temperature")
    plt.ylabel("Heat Capacity")
    plt.errorbar(temp, heat_cap_lst, yerr=std)
    plt.savefig('heat_capacity_g.png')
    plt.show()

       
#Kawasaki dynamics         
if dynamic == "K":
    #set to ground state configuration (one half all up, the other all down)
    lattice = np.ones((Nx,Ny), dtype=int)
    lattice[int(N/2):int(N), 0:int(N)] = -1
    print(lattice)
    #initial temperature
    T = Ti
    while T <= Tf:
        #append temperature to list
        temp.append(T)
        print(T)
        #initiate energy lists
        energy_all = []   
        for n in range(steps):
            for a in range(Nx*Ny):
                
                #randomly generate two sites
                x1 = random.randint(0, Nx-1)
                y1 = random.randint(0, Ny-1)
                x2 = random.randint(0, Nx-1)
                y2 = random.randint(0, Ny-1)
                
                #determine new spins to assign to random coordinates after swap
                new_spin1 = lattice[x2,y2]
                new_spin2 = lattice[x1,y1]
                
                #if the spins are the same no change in energy
                if new_spin1 == new_spin2:
                    delta_energy = 0
                #if spins are different energy changes
                else:
                    delta_energy = -2*energy(lattice,x1,y1)-2*energy(lattice,x2,y2)
                    #nearest neighbours energy caluclation 
                    #Kawasaki energy calculation (Glauber at 2 locations with correction) 
                    if x1+1 == x2 and y1 == y2: 
                        delta_energy += 4
                    if x1-1 == x2 and y1 == y2: 
                        delta_energy += 4
                    if x1 == x2 and y1+1 == y2: 
                        delta_energy += 4
                    if x1 == x2 and y1-1 == y2:
                        delta_energy += 4
                    #non-nearest neighbour energy calculation     
                    else:
                        #treat as Glauber at 2 locations
                        delta_energy = delta_energy
                
                #metrapolis algorithm - determine when to allow spin reverse
                if delta_energy < 0:
                    lattice[x1,y1] = new_spin1
                    lattice[x2,y2] = new_spin2
                else:
                    prob = math.exp(-delta_energy/T)
                    if random.random() <= prob:
                        lattice[x1,y1] = new_spin1
                        lattice[x2,y2] = new_spin2
                        
            #calculate energy every 10 steps (magnetisation is constant)    
            if n >= 100 and n % 10 == 0:
                energy_tot = 0
                for i in range(Nx):
                    for j in range(Ny):
                        energy_tot += energy(lattice,i,j)
                energy_all.append(0.5*energy_tot)
            
            if animate == 'Y':
                if(n % 100 == 0):        
                    f = open('spins_k.dat','w')
                    for i in range(Nx):
                        for j in range(Ny):
                            f.write('%d %d %lf\n'%(i,j,lattice[i,j]))
                    f.close()
                    #show animation
                    plt.cla()
                    im=plt.imshow(lattice, animated=True)
                    plt.draw()
                    plt.pause(0.0001)
        
        #calculate average energy and heat capacity across temperatures
        energy_all_sq = np.array(energy_all)**2
        avg_energy = np.mean(energy_all)
        avg_energy_lst.append(avg_energy)
        heat_cap = (np.mean(energy_all_sq)-(avg_energy**2))/((N**2)*(T**2))
        heat_cap_lst.append(heat_cap)
        
        #bootstrap error
        m = 0
        heat_cap_s_lst = []
        while m <= 1000:
            #get sample (with replacment)
            s = np.random.choice(energy_all, len(energy_all))
            avg_s = np.mean(s)
            avg_s_sq = np.array(s)**2
            #calculate heat capacity for sample
            heat_cap_s = (np.mean(avg_s_sq)-(avg_s**2))/((N**2)*(T**2))
            #append heat capacity to a list of all sample heat capacities
            heat_cap_s_lst.append(heat_cap_s)
            m += 1
            #calculate variance and standard deviation    
        var = (np.mean(np.array(heat_cap_s_lst)**2)-np.mean(heat_cap_s_lst)**2)
        std = var**(1/2) 
        errors.append(std)
         
        #increase temperature
        T += t
        

    
    #write data to files
    f = open('energy_k.dat','w')
    for v in range(len(temp)):
            f.write(str(temp[v]) + ": " + str(avg_energy_lst[v]) + "\n")
    f.close()
    
    f = open('heat_capacity_k.dat','w')
    for v in range(len(temp)):
            f.write(str(temp[v]) + ": " + str(heat_cap_lst[v]) + " +/- " + str(errors[v]) + "\n")
    f.close()       
                    
    #produce plots    
    plt.title("Energy change with temperature")
    plt.xlabel("Temperature")
    plt.ylabel("Energy")
    plt.plot(temp, avg_energy_lst)
    plt.savefig('energy_k.png')
    plt.show()
    
    
    plt.title("Heat Capacity change with temperature")
    plt.xlabel("Temperature")
    plt.ylabel("Heat Capacity")
    plt.errorbar(temp, heat_cap_lst, yerr=errors)
    plt.savefig('heat_capacity_k.png')
    plt.show()

       
#show updated lattice
#fig = plt.figure()
#im = plt.imshow(lattice, animated = True)


