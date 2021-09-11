#This will just be some coding to get used to working with python again
#Computing aims: 
# 1: Modular
# 2: Structures / Classes
# 3: Reading and Writing to files
# 4: Drawing Graphs
# 5: Creating GIFS

#Basic principle of this program is to create functions for calculating characteristics of the particles in another file
#called functions, create a CSV file of all the data as going along then recalling a function called plot, in the plot file to create 
#the gif and graph
#Imports

from matplotlib import pyplot as plt
from matplotlib import animation
import numpy as np
import Functions # Module File:)
import random
import Plot
import os
import math

#Setting up variable to hold information
particles = []

#(self, mass, radius, position, velocity, acceleration, force, KE, PE)
#Creating the different particles in the model 

#This is to randomly choose the sign on positions so i don't have to go over 0
plusminus = [-1,1]
for i in range (Functions.particle_no -1):
 
    pos = random.uniform(6371000, 8000000)
    
    Xx = random.uniform(-7000000,7000000) * random.choice(plusminus)
    Xy_range = ((pos**2 - Xx**2)**(1/2)).real
    Xy = random.uniform(-Xy_range,Xy_range)
    Xz = ((pos**2 - Xx**2 - Xy**2)**(1/2) * random.choice(plusminus)).real
        
            
    velo = random.uniform(7000,8000)
        
    Vx = random.uniform(-8000,8000)
    Vy_range = ((velo**2 - Vx**2)**(1/2)).real
    Vy = random.uniform(-Vy_range,Vy_range)
    Vz = ((velo**2 - Vx**2 - Vy**2)**(1/2)  * random.choice(plusminus)).real
       
       
    M  = random.uniform(1*10**10,1*10**5)
    r  = random.uniform(1,10)
    particles.append(Functions.Particle(M,r, np.array([Xx,Xy,Xz]), np.array([Vx,Vy,Vz]), np.array([0,0,0]), 0, 0, 0))

#Saturn
particles.append(Functions.Particle(5.683*10**24,50 ,np.array([0,0,0]), np.array([0,0,0]), np.array([0,0,0]), 0, 0, 0)) #This is the fat planet


#The main module
def main():
    #Opening the file to write to
    file = open("DataDump.csv", "a")
    

    #Grabbing inputs for number of iterations and iterations per frame
    n = Functions.Integer_Input_Positive("How many iterations? ")
    while True:

        Plot.iterations_per_frame = Functions.Integer_Input_Positive("How many iterations per frame? ")
        if Plot.iterations_per_frame > n:
            print("Enter a value less than the number of iterations ")
            continue
        else:
            break
        
    i=0  #Step counter
    t = 0 #Track the time     

    #Initialising Values
    Functions.KineticEnergy(particles)
    Functions.Gravity(particles)

    #FileWritingTime!

    #File layout
    # Step, Time, X_x, X_y, X_z, V_x, V_y, V_z, A_x, A_y, A_z, KE, PE
    file.write("Step, Time, X_x, X_y, X_z, V_x, V_y, V_z, A_x, A_y, A_z, KE, PE, radius") #So I know the layout when I read the file in the future
    

    while i <= n:
        
        file.write("\n") #New line for new time
        initial = ( str(i) +"," + str(t) + ",")
        file.write(str(initial)) #Writing File and step to the file
        
        for particle in particles:
                       
            #Writing Displacement
            X = str((particle.position[0])) + "," + str((particle.position[1])) + ","+ str((particle.position[2])) + ","

            file.write(X)
        
            #Writing Velocity 
            V = str((particle.velocity[0])) + "," + str((particle.velocity[1])) + "," + str((particle.velocity[2])) + ","
            file.write(V)

            #Writing Acceleration
            A = str((particle.acceleration[0])) + "," + str((particle.acceleration[1])) + "," + str((particle.acceleration[2])) + ","
            file.write(A)
            
            #Writing KE and PE
            KE_PE = str(particle.KE) + "," + str(particle.PE) + ","
            file.write(KE_PE)

            #Writing radius
            radius = str(particle.radius) + ","
            file.write(radius)
            
        #Updates position and velocity
        Functions.Step(particles)

        #Updates the force,acceleration, KE,PE
        Functions.Gravity(particles)
        Functions.KineticEnergy(particles)
        
        #Counting functions
        t += Functions.dt
        i+=1

    file.close()
    return(None)

#Running the functions
main()
Plot.PlottingAll()
Plot.Plotting3D()
Plot.Plotting2D()
Plot.PlottingEnergy()
#Removes the csv file to leave only the gif behind 
#Take out when i've finished the testing phase 
os.remove("DataDump.csv")
