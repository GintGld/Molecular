from Constants import *
from Classes import *
from Functions import *
import scipy as sp
import matplotlib.pyplot as plt
import numpy as np
import random as rn
# Creating constants
Create_new_patricles = True
# Modeling constants
Do_modeling = True
Save_modeling = True
Record_energy = True
Record_velocity = True
Write_VMD = True
# Analysis constants
Graph_velocity = False
Energy_distribution = False
Graph_energy = True

def main():
    if Create_new_patricles:
        '''
                Создание массива частиц
        '''
        Generate_particles()
        print('Particles were generated.')
    if Do_modeling:
        '''
                Моделирование движения
        '''
        if not(Number_of_particles >= Particles_in_one_row and Number_of_particles <= Particles_in_one_row ** 3):
            print('Wrong constants. Fix Number_of_particles and Particles_in_one_row.')
            return 0
        input_P = open('initial_conditions.txt', 'r')
        Particles = []
        for i in range(Number_of_particles):
            s = input_P.readline()
            s = s.rstrip()
            s = s.split(',')
            Particles.append(particle(float(s[0]), float(s[1]), float(s[2]), 
            float(s[3]), float(s[4]), float(s[5])))
        input_P.close()
        if Record_energy:
            output_E = open('Energy.txt', 'w')
        if Write_VMD:
            output_VMD = open('Coordinates_VMD.txt', 'w')
        V_mean = [0] * len(Particles)
        for i in range(1, Number_of_steps + 1):
            if i + 7 >= Number_of_steps:
                V = Get_velocity(Particles)
                for i in range(0, len(V_mean)):
                    V_mean[i] += V[i] / 7
            if i == 1:
                print('Modeling started...')
            if Record_energy:
                print(Get_energy(Particles, True), file=output_E)
            Update_positions(Particles)
            if i % Step_of_count == 0:
                print(str(i) + ' of ' + str(Number_of_steps) + ' steps were modeled...')
            if Write_VMD:
                print(Number_of_particles, file=output_VMD)
                print('', file=output_VMD)
                for part in Particles:
                    print('1 ' + str(part.x) + ' ' + str(part.y) + ' ' + str(part.z), file=output_VMD)
        if Record_energy:
            print(Get_energy(Particles, True), file=output_E)
            output_E.close()
        if Write_VMD:
            output_VMD.close()
            print('Coordinates for VMD were recorded.')
        if Record_velocity:
            output_V = open('Velocity.txt', 'w')
            V = Get_velocity(Particles)
            for i in range(len(V)):
                print(V[i], file=output_V)
            output_V.close()
        print('Modeling ended.')
        if Save_modeling:
            output_cond = open('initial_conditions.txt', 'w')
            for part in Particles:
                print(str(part.x) + ',' + str(part.y) + ',' + str(part.z) + ',' + 
                str(part.vx) + ',' + str(part.vy) + ',' + str(part.vz), file=output_cond)
            output_cond.close()
            print('Final conditions was recorded.')
        if Record_energy:
            print('Energy was recorded.')
        if Record_velocity:
            print('Velocity was recorded.')
    if Graph_velocity:
        '''
                График распределения скоростей
        '''
        input_V = open('velocity.txt', 'r')
        V = []
        for i in range(Number_of_particles):
            s = input_V.readline()
            s = s.rstrip()
            V.append(float(s))
        input_V.close()
        if (len(V) != Number_of_particles):
            print('Wrong information, do modeling.')
        #else:
            #plt.plot(V)
        r = 2 * (np.percentile(V_mean, 75) - np.percentile(V_mean, 25)) / (Number_of_particles ** (1 / 3))
        plt.hist(V_mean, bins=np.arange(0, max(V_mean) + r, r))
        plt.show()
    if Energy_distribution:
        '''
                Проверка нормального распределения скоростей
        '''
        input_V = open('Energy.txt', 'r')
        V = []
        for i in range(Number_of_steps):
            s = input_V.readline()
            s = s.rstrip()
            V.append(float(s))
        input_V.close()
        if len(V) != Number_of_steps:
            print('Wrong information, do modeling.')
        else:
            Dispersion(V)
    if Graph_energy:
        input_E = open('Energy.txt', 'r')
        E = []
        x = []
        for i in range(Number_of_steps):
            s = input_E.readline()
            s = s.strip()
            E.append(float(s))
            x.append(i + 1)
        input_E.close()
        plt.plot(x, E)
        plt.show()
    

if __name__ == "__main__":
    main()