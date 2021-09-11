#This is the module where I define all the functions for the actual program... yay!

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation

particle_no = 50
G = 6.67430*(10**-11)
dt = 1

#Creating the class for each particle/object
class Particle:
    def __init__(self, mass, radius, position, velocity, acceleration, force, KE, PE): #Setting up the class for the general particle (all in standard SI units)
        self.mass         = mass
        self.radius       = radius
        self.position     = position
        self.velocity     = velocity
        self.acceleration = acceleration
        self.force        = force
        self.KE           = KE
        self.PE           = PE

#Calculating the gravitational force and potential on a particle/object

#I want nested loops to calculate the force between each particle if I want to expand this to N particles
#The general idea of this algorithm is to calcualte the force vectors / potential from all particles acting on a single particle
# add them up and then continue with that for each particle
#Returns arrays of the forces on all the particles and then the potential energies

def Gravity(particles): 
    
    for particle1 in particles:
        
        total_force = np.array([0,0,0], float)   #The variables added to, to find the totals
        total_potential = 0

        for particle2 in particles:
            if particle1 == particle2: #Just so we don't calcualte the force/potential from gravity of a particle on itself 
                continue

            r         = abs(np.linalg.norm(particle1.position -particle2.position))
            force     = (G*particle1.mass*particle2.mass)/r**2
            potential = -(G*particle1.mass*particle2.mass)/r
            direction = (particle2.position-particle1.position)/r
            
            total_potential += potential
            total_force     += direction*force
    
        particle1.force = total_force
        particle1.acceleration = total_force / particle1.mass
        particle1.PE = total_potential

    return None

#Kinetic Energy of the particles (It just works them all out)
def KineticEnergy(particles):
    for particle in particles:
        velocity_mag = np.linalg.norm(particle.velocity)
        KE = 0.5 * particle.mass * ((velocity_mag)**2)
        particle.KE = KE
    return



#Step of the position and velocity of the particle/object
def Step(particles):
    for particle in particles:
        initial_pos = particle.position
        initial_velo = particle.velocity
        
        particle.position = initial_pos + particle.velocity*dt + (1/2)*particle.acceleration*(dt**2)
        particle.velocity = initial_velo + particle.acceleration*dt
        #k_1 = particle.velocity
        #k_2 = particle.velocity*(dt/2) + particle.velocity*dt*k_1/2 
        #k_3 = particle.velocity*(dt/2) + particle.velocity*dt*k_2/2
        #k_4 = particle.velocity*dt + particle.velocity*dt*k_3 
        #particle.position = initial_pos + (dt/6)*(k_1 +2*k_2 + 2*k_3 + k_4)

        #k_1 = particle.acceleration
        #k_2 = particle.acceleration + particle.acceleration*dt*k_1/2 
        #k_3 = particle.acceleration + particle.acceleration*dt*k_2/2
        #k_4 = particle.acceleration + particle.acceleration*dt*k_3 
        #particle.velocity = initial_velo + (dt/6)*(k_1 +2*k_2 + 2*k_3 + k_4)


    return None

def Integer_Input_Positive(message):
    while True:
        try:
            userInput = int(input(message))       
        except ValueError:
            print("Not an integer! Try again.")
            continue
        else:
            if userInput <= 0:
                print("Enter a value greater than 0.")
            
            else:    
                return userInput 
                break 
