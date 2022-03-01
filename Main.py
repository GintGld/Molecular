from Constants import *
from Classes import *
from Functions import *
import matplotlib.pyplot as plt
import numpy as np
import random as rn

Particles = []

def Create():
    '''
            Создание массива частиц
    '''
    Generate_particles()
    print('Particles were generated.')

def Modeling():
    '''
            Моделирование движения
    '''
    if not(Number_of_particles >= Particles_in_one_row and Number_of_particles <= Particles_in_one_row ** 3):
        print('Wrong constants. Fix Number_of_particles and Particles_in_one_row.')
        return 0
    input_P = open('initial_conditions.txt', 'r')
    for i in range(Number_of_particles):
        s = input_P.readline()
        s = s.rstrip()
        s = s.split(',')
        Particles.append(particle(float(s[0]), float(s[1]), float(s[2]), 
        float(s[3]), float(s[4]), float(s[5])))
    input_P.close()
    output_E = open('Energy.txt', 'w')
    output_K = open('Kinetic_energy.txt', 'w')
    output_P = open('Potential_energy.txt', 'w')
    output_M = open('Momentum.txt', 'w')
    output_VMD = open('Coordinates_VMD.txt', 'w')
    V_mean = [0] * len(Particles)
    print('Modeling started...')

    for i in range(1, Number_of_steps + 1):
        if i + 100 >= Number_of_steps:
            V = Get_velocity(Particles)
            for i in range(0, len(V_mean)):
                V_mean[i] += V[i] / 100
        print(Get_energy(Particles, True), file=output_E)
        print(K_energy(Particles), file=output_K)
        print(P_energy(Particles), file=output_P)
        P = Get_momentum(Particles)
        print(str(P[0]) + ',' + str(P[1]) + ',' + str(P[2]), file=output_M)
        Update_positions(Particles)
        print(Number_of_particles, file=output_VMD)
        print('', file=output_VMD)
        for part in Particles:
            print('1 ' + str(part.x) + ' ' + str(part.y) + ' ' + str(part.z), file=output_VMD)
        if i % Step_of_count == 0:
            print(str(i) + ' of ' + str(Number_of_steps) + ' steps were modeled...')
            
    print(Get_energy(Particles, True), file=output_E)
    print(K_energy(Particles), file=output_K)
    print(P_energy(Particles), file=output_P)
    output_E.close()
    output_K.close()
    output_P.close()
    output_VMD.close()

    P = Get_momentum(Particles)
    print(str(P[0]) + ',' + str(P[1]) + ',' + str(P[2]), file=output_M)
    output_M.close()

    output_V = open('Velocity.txt', 'w')
    for i in range(len(V_mean)):
        print(V_mean[i], file=output_V)
    output_V.close()

    print(Number_of_steps + ' of ' + Number_of_steps + ' were modeled...')
    print('Modeling ended.')

def Save_Modeling():
    output_cond = open('initial_conditions.txt', 'w')
    for part in Particles:
        print(str(part.x) + ',' + str(part.y) + ',' + str(part.z) + ',' + 
        str(part.vx) + ',' + str(part.vy) + ',' + str(part.vz), file=output_cond)
    output_cond.close()

def Velocity():
    '''
            График распределения скоростей
    '''
    input_V = open('velocity.txt', 'r')
    V = []
    for i in range(Number_of_particles):
        s = input_V.readline()
        s = s.rstrip()
        if (s == ''): break
        V.append(float(s))
    input_V.close()
    if (len(V) != Number_of_particles):
        print('Wrong information, do modeling.')
    r = 2 * (np.percentile(V, 75) - np.percentile(V, 25)) / (Number_of_particles ** (1 / 3))
    r = (np.percentile(V, 75) - np.percentile(V, 25)) / (Number_of_particles ** (1 / 3))
    plt.hist(V, bins=np.arange(0, max(V) + r, r))
    plt.show()

def Energy():
    input_E = open('Energy.txt', 'r')
    input_EK = open('Kinetic_energy.txt', 'r')
    input_EP = open('Potential_energy.txt', 'r')
    E = []
    E_K = []
    E_P = []
    for i in range(Number_of_steps):
        s = input_E.readline()
        s = s.strip()
        if s == '': break
        E.append(float(s))
        s = input_EK.readline()
        s = s.strip()
        if s == '': break
        E_K.append(float(s))
        s = input_EP.readline()
        s = s.strip()
        if s == '': break
        E_P.append(float(s))
    input_E.close()
    input_EK.close()
    input_EP.close()
    x1 = np.arange(0, len(E), 1)
    x2 = np.arange(0, len(E_K), 1)
    x3 = np.arange(0, len(E_P), 1)
    plt.plot(x1, E, label='Full Energy')
    plt.plot(x2, E_K, label='Kinetic Energy')
    plt.plot(x3, E_P, label='Potential Energy')
    plt.legend(loc='best', fontsize=12)
    plt.xlabel('Number of step', fontsize=14)
    plt.ylabel('Energy', fontsize=14)
    plt.grid(True)
    plt.show()

def Momentum():
    input_M = open('Momentum.txt', 'r')
    P = []
    for i in range(Number_of_steps):
        s = input_M.readline()
        s = s.strip()
        s = s.split(',')
        if s == '': break
        P.append(np.sqrt(float(s[0]) ** 2 + float(s[1]) ** 2 + float(s[2]) ** 2))
    x = np.arange(len(P))
    plt.plot(x, P, label='Full')
    plt.legend(loc='best', fontsize=12)
    plt.xlabel('Number of step', fontsize=14)
    plt.ylabel('Momentum', fontsize=14)
    plt.show()

if __name__ == "__main__":
    Create()
    Modeling()
    Save_Modeling()
    Velocity()
    Energy()
    #Momentum()