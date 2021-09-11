#File to plot all the data

import matplotlib.pyplot as plt
import numpy as np
import Functions
import csv
import imageio
import os
from pygifsicle import optimize
from mpl_toolkits import mplot3d
import PIL
import math

iterations_per_frame = 0

def PlottingAll():
    #Setting up variables to create the GIF later
    gif_name = 'All' 
    n_frames=10
    marker_size = 25

    #Setting up the file to be read
    file = open("DataDump.csv", "r")
    reader = csv.reader(file, delimiter = ',', skipinitialspace=True)

    #Reads all the rows of the files, and how many there are
    rows = list(reader)
    row_numb = int(len(rows)) - 1

    filenames = []

    #Setting up lists to contain energies of each row of data
    KE = []
    PE = []
    E  = []
    Time = []
    for i in range (row_numb):
        
        #Limit for 3D and 2D graphs
        sq = 10000000

        n=0 #Reset counter
        fig = plt.figure(figsize=(20,20))
        
        #Creating 3 Subplots, 1 3D, another 2D, another to scatter the energy
        ax1 = plt.subplot(2,2,1,projection='3d')
        ax1.set_title("3D XYZ")

        #Setting Limits on 3D Graph
        ax1.set_xlim([-1*sq,sq])
        ax1.set_ylim([-1*sq,sq])
        ax1.set_zlim([-1*sq,sq])

        #Labelling Axis 3D Graph
        ax1.set_xlabel('X Label')
        ax1.set_ylabel('Y Label')
        ax1.set_zlabel('Z Label')

        #2D XY Plane Graph
        ax2 = plt.subplot(2,2,2)
        ax2.set_title("2D XY Plane")
        ax2.set_aspect('equal', adjustable='box')

        #limits
        ax2.set_xlim([-1*sq,sq])
        ax2.set_ylim([-1*sq,sq])

        #Labels
        ax2.set_xlabel('X Label')
        ax2.set_ylabel('Y Label')

        #Energy Graph
        ax3 = plt.subplot(2,1,2)
        ax3.set_title("Total Energy, KE, PE against T")

        #Limits
        total_time = float(((rows[row_numb])[1]))
        
        
        #Going to find the max E and set 3*Max E as the axis plusminus limits for this plot
        E_max =0

        for a in range (row_numb):
            KE_temp = 0
            PE_temp = 0 
            while n < Functions.particle_no:

                KE_temp += float(((rows[1+a])[11+12*n]))
                PE_temp += float(((rows[1+a])[12+12*n]))
                n+=1

            E_temp = KE_temp + PE_temp
            if abs(E_temp) > abs(E_max):
                E_max = E_temp
            
            if abs(KE_temp) > abs(E_max):
                E_max = KE_temp

            if abs(PE_temp) > abs(E_max):
                E_max = PE_temp
        
        ax3.set_xlim([0,total_time])
        ax3.set_ylim([-1*abs(E_max)*3,abs(E_max)*3])

        KE_temp =0
        PE_temp =0
        E_temp =0
        n=0
        filenames.append(str(i) + '.png')
        while n < Functions.particle_no:
            
            f = iterations_per_frame*i  #After how many rows does it update
            
            #Reading the data from each particle 
            
            Xx_temp = float(((rows[1+f])[2+12*n]))
            Xy_temp = float(((rows[1+f])[3+12*n]))
            Xz_temp = float(((rows[1+f])[4+12*n]))

            KE_temp += float(((rows[1+f])[11+12*n]))
            PE_temp += float(((rows[1+f])[12+12*n]))
            
            radius_temp = float((rows[1+i])[13+12*n])


            #Plotting the points
            ax1.scatter(Xx_temp,Xy_temp,Xz_temp,s=radius_temp)
            ax2.scatter(Xx_temp,Xy_temp,s=radius_temp)
            
            #Plotting directions of velocity and acceleration
            #ax.quiver(Xx_temp, Xy_temp, Xz_temp, Vx_temp ,Vy_temp , Vz_temp)
            #ax.quiver(Xx_temp, Xy_temp, Xz_temp, Ax_temp ,Ay_temp , Az_temp)
        
            n+=1
        current_time = float(((rows[1+f])[1]))
        E_temp = abs(KE_temp + (PE_temp)/2)
        KE.append(KE_temp)
        PE.append(PE_temp)
        E.append(E_temp)
        Time.append(current_time)
        ax3.plot(Time,KE, label = "KE")
        ax3.plot(Time,PE, label = "PE")
        ax3.plot(Time,E, label = "Total E")
        ax3.legend()

        plt.savefig(str(i) + '.png', dpi=80)
        plt.close()
        if f == (row_numb - iterations_per_frame) or f > (row_numb - iterations_per_frame):
                break
        
    #Building the gif
    with imageio.get_writer(f'{gif_name}.gif', mode='I') as writer:
        for filename in filenames:
            image = imageio.imread(filename)
            writer.append_data(image)

    #Remove all the created images after the GIF has been made
    for filename in set(filenames):
       os.remove(filename)

    #https://towardsdatascience.com/basics-of-gifs-with-pythons-matplotlib-54dd544b6f30

    file.close()

def Plotting3D():
    #Some gif stuff
    gif_name = '3D' 
    n_frames=10
    marker_size = 25

    #Setting up the file to be read
    file = open("DataDump.csv", "r")
    reader = csv.reader(file, delimiter = ',', skipinitialspace=True)

    #Reads all the rows of the files, and how many there are
    rows = list(reader)
    row_numb = int(len(rows)) - 1

    filenames = []

    #Setting up lists to contain energies of each row of data
    KE = []
    PE = []
    E  = []
    Time = []
    for i in range (row_numb):
        
        #Limit for 3D and 2D graphs
        sq = 10000000

        n=0 #Reset counter
        fig = plt.figure(figsize=(10,10))
        
        #Creating 3 Subplots, 1 3D, another 2D, another to scatter the energy
        ax1 = plt.subplot(projection='3d')
        ax1.set_title("3D XYZ")

        #Setting Limits on 3D Graph
        ax1.set_xlim([-1*sq,sq])
        ax1.set_ylim([-1*sq,sq])
        ax1.set_zlim([-1*sq,sq])

        #Labelling Axis 3D Graph
        ax1.set_xlabel('X Label')
        ax1.set_ylabel('Y Label')
        ax1.set_zlabel('Z Label')

        
        KE_temp =0
        PE_temp =0
        E_temp =0
        n=0
        filenames.append(str(i) + '.png')
        while n < Functions.particle_no:
            
            f = iterations_per_frame*i  #After how many rows does it update
            
            #Reading the data from each particle 
            
            Xx_temp = float(((rows[1+f])[2+12*n]))
            Xy_temp = float(((rows[1+f])[3+12*n]))
            Xz_temp = float(((rows[1+f])[4+12*n]))
            
            KE_temp += float(((rows[1+f])[11+12*n]))
            PE_temp += float(((rows[1+f])[12+12*n]))
            

            
            radius_temp = float((rows[1+i])[13+12*n])


            #Plotting the points
            ax1.scatter(Xx_temp,Xy_temp,Xz_temp,s=radius_temp)

            #Plotting directions of velocity and acceleration
            #ax.quiver(Xx_temp, Xy_temp, Xz_temp, Vx_temp ,Vy_temp , Vz_temp)
            #ax.quiver(Xx_temp, Xy_temp, Xz_temp, Ax_temp ,Ay_temp , Az_temp)
        
            n+=1
        current_time = float(((rows[1+f])[1]))
        E_temp = abs(KE_temp + (PE_temp)/2)
        KE.append(KE_temp)
        PE.append(PE_temp)
        E.append(E_temp)
        Time.append(current_time)
        

        plt.savefig(str(i) + '.png', dpi=80)
        plt.close()
        if f == (row_numb - iterations_per_frame) or f > (row_numb - iterations_per_frame):
                break
        
    #Building the gif
    with imageio.get_writer(f'{gif_name}.gif', mode='I') as writer:
        for filename in filenames:
            image = imageio.imread(filename)
            writer.append_data(image)

    #Remove
    for filename in set(filenames):
       os.remove(filename)

    #https://towardsdatascience.com/basics-of-gifs-with-pythons-matplotlib-54dd544b6f30

    file.close()

def Plotting2D():
    #Some gif stuff
    gif_name = '2D' 
    n_frames=10
    marker_size = 25

    #Setting up the file to be read
    file = open("DataDump.csv", "r")
    reader = csv.reader(file, delimiter = ',', skipinitialspace=True)

    #Reads all the rows of the files, and how many there are
    rows = list(reader)
    row_numb = int(len(rows)) - 1

    filenames = []

    #Setting up lists to contain energies of each row of data
    KE = []
    PE = []
    E  = []
    Time = []
    for i in range (row_numb):
        
        #Limit for 3D and 2D graphs
        sq = 10000000

        n=0 #Reset counter
        fig = plt.figure(figsize=(10,10))
        
        #Creating 3 Subplots, 1 3D, another 2D, another to scatter the energy  

        #2D XY Plane Graph
        ax2 = plt.subplot()
        ax2.set_title("2D XY Plane")
        ax2.set_aspect('equal', adjustable='box')

        #limits
        ax2.set_xlim([-1*sq,sq])
        ax2.set_ylim([-1*sq,sq])

        #Labels
        ax2.set_xlabel('X Label')
        ax2.set_ylabel('Y Label')
        
        KE_temp =0
        PE_temp =0
        E_temp =0
        n=0
        filenames.append(str(i) + '.png')
        while n < Functions.particle_no:
            
            f = iterations_per_frame*i  #After how many rows does it update
            
            #Reading the data from each particle 
            
            Xx_temp = float(((rows[1+f])[2+12*n]))
            Xy_temp = float(((rows[1+f])[3+12*n]))

            KE_temp += float(((rows[1+f])[11+12*n]))
            PE_temp += float(((rows[1+f])[12+12*n]))
                        
            radius_temp = float((rows[1+i])[13+12*n])

            #Plotting the points
            ax2.scatter(Xx_temp,Xy_temp,s=radius_temp)
            
            #Plotting directions of velocity and acceleration
            #ax.quiver(Xx_temp, Xy_temp, Xz_temp, Vx_temp ,Vy_temp , Vz_temp)
            #ax.quiver(Xx_temp, Xy_temp, Xz_temp, Ax_temp ,Ay_temp , Az_temp)
        
            n+=1
        current_time = float(((rows[1+f])[1]))
        E_temp = abs(KE_temp + (PE_temp)/2)
        KE.append(KE_temp)
        PE.append(PE_temp)
        E.append(E_temp)
        Time.append(current_time)
        
        plt.savefig(str(i) + '.png', dpi=80)
        plt.close()
        if f == (row_numb - iterations_per_frame) or f > (row_numb - iterations_per_frame):
                break
        
    #Building the gif
    with imageio.get_writer(f'{gif_name}.gif', mode='I') as writer:
        for filename in filenames:
            image = imageio.imread(filename)
            writer.append_data(image)

    #Remove
    for filename in set(filenames):
        os.remove(filename)

    #https://towardsdatascience.com/basics-of-gifs-with-pythons-matplotlib-54dd544b6f30

    file.close()

def PlottingEnergy():
    #Some gif stuff
    gif_name = 'Energy' 
    n_frames=10
    marker_size = 25

    #Setting up the file to be read
    file = open("DataDump.csv", "r")
    reader = csv.reader(file, delimiter = ',', skipinitialspace=True)

    #Reads all the rows of the files, and how many there are
    rows = list(reader)
    row_numb = int(len(rows)) - 1

    filenames = []

    #Setting up lists to contain energies of each row of data
    KE = []
    PE = []
    E  = []
    Time = []
    for i in range (row_numb):
        
        #Limit for 3D and 2D graphs
        sq = 10000000

        n=0 #Reset counter
        fig = plt.figure(figsize=(10,10))
        
        #Creating 3 Subplots, 1 3D, another 2D, another to scatter the energy
    
        #Energy Graph
        ax3 = plt.subplot()
        ax3.set_title("Total Energy, KE, PE against T")

        #Limits
        total_time = float(((rows[row_numb])[1]))
        
        #Going to find the max E and set 1.1*Max E as the Limit for this plot ITS FOR EACH ROW NOT PARTICLE, i,e all particles
        E_max =0

        for a in range (row_numb):
            KE_temp = 0
            PE_temp = 0 
            while n < Functions.particle_no:

                KE_temp += float(((rows[1+a])[11+12*n]))
                PE_temp += float(((rows[1+a])[12+12*n]))
                n+=1

            E_temp = KE_temp + PE_temp
            if abs(E_temp) > abs(E_max):
                E_max = E_temp
            
            if abs(KE_temp) > abs(E_max):
                E_max = KE_temp

            if abs(PE_temp) > abs(E_max):
                E_max = PE_temp
        
        ax3.set_xlim([0,total_time])
        ax3.set_ylim([-1*abs(E_max)*3,abs(E_max)*3])

        
        KE_temp =0
        PE_temp =0
        E_temp =0
        n=0
        filenames.append(str(i) + '.png')
        while n < Functions.particle_no:
            
            f = iterations_per_frame*i  #After how many rows does it update
            
            #Reading the data from each particle 
            #         
            KE_temp += float(((rows[1+f])[11+12*n]))
            PE_temp += float(((rows[1+f])[12+12*n]))
            
            #Plotting the points        
            #Plotting directions of velocity and acceleration
            #ax.quiver(Xx_temp, Xy_temp, Xz_temp, Vx_temp ,Vy_temp , Vz_temp)
            #ax.quiver(Xx_temp, Xy_temp, Xz_temp, Ax_temp ,Ay_temp , Az_temp)
        
            n+=1
        current_time = float(((rows[1+f])[1]))
        E_temp = abs(KE_temp + (PE_temp)/2)
        KE.append(KE_temp)
        PE.append(PE_temp)
        E.append(E_temp)
        Time.append(current_time)
        ax3.plot(Time,KE, label = "KE")
        ax3.plot(Time,PE, label = "PE")
        ax3.plot(Time,E, label = "Total E")
        ax3.legend()

        plt.savefig(str(i) + '.png', dpi=80)
        plt.close()
        if f == (row_numb - iterations_per_frame) or f > (row_numb - iterations_per_frame):
                break
        
    #Building the gif
    with imageio.get_writer(f'{gif_name}.gif', mode='I') as writer:
        for filename in filenames:
            image = imageio.imread(filename)
            writer.append_data(image)

    #Remove
    for filename in set(filenames):
        os.remove(filename)

    #https://towardsdatascience.com/basics-of-gifs-with-pythons-matplotlib-54dd544b6f30

    file.close()
