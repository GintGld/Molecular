from Constants import *
from Classes import *
import numpy as np
import random as rn

def Generate_particles():
    '''
        Генерирует массив частиц, расположенных на расстоянии не больших четверти шага сетки 
        с рандомными  скоростями, ограниченные константой  V_MAX
        Возвращает массив с элементами класса particle
    '''
    vx = np.random.normal(0, 0.5, Number_of_particles)
    vy = np.random.normal(0, 0.5, Number_of_particles)
    vz = np.random.normal(0, 0.5, Number_of_particles)
    output = open('initial_conditions.txt', 'w')
    Particles = []
    counter = 0
    for i in range(1, Particles_in_one_row + 1):
        for j in range(1, Particles_in_one_row + 1):
            for k in range(1, Particles_in_one_row + 1):
                if (len(Particles) == Number_of_particles):
                    break
                Particles.append(particle(
                (i + (rn.random() - 0.5) * 2 * Rad_of_rand_gen) * Size_of_cell, 
                (j + (rn.random() - 0.5) * 2 * Rad_of_rand_gen) * Size_of_cell, 
                (k + (rn.random() - 0.5) * 2 * Rad_of_rand_gen) * Size_of_cell,
                vx[counter], vy[counter], vz[counter]))
    for part in Particles:
        print(str(part.x) + ',' +  str(part.y) + ',' +  str(part.z) + ',' +  
        str(part.vx) + ',' +  str(part.vy) + ',' +  str(part.vz),  file=output)
    output.close()

def Update_acceleration(Particles):
    '''
        Производит покоординатный расчет ускорений для каждой частицы
    '''
    for part in Particles:
        part.prev_ax = part.ax
        part.prev_ay = part.ay
        part.prev_az = part.az
        part.ax = 0
        part.ay = 0
        part.az = 0
        for obj in Particles:
            if part != obj:
                x, y, z = part.nearest_reflection(obj)
                l = Distance(part, x, y, z)
                if l >= 2.5: continue
                part.ax += 24 * (part.x - x) * (2 * l ** (-14) - l ** (-8))
                part.ay += 24 * (part.y - y) * (2 * l ** (-14) - l ** (-8))
                part.az += 24 * (part.z - z) * (2 * l ** (-14) - l ** (-8))
    
def Update_velocity(Particles):
    '''
        Рассчет моментальной скорости каждой частицы
    '''
    for part in Particles:
        if part.prev_ax == 0 and part.prev_ay == 0 and part.prev_az == 0:
            part.prev_ax = part.ax
            part.prev_ay = part.ay
            part.prev_az = part.az
        part.vx += 0.5 * (part.prev_ax + part.ax) * dt
        part.vy += 0.5 * (part.prev_ay + part.ay) * dt
        part.vz += 0.5 * (part.prev_az + part.az) * dt

def Get_energy(Particles, type):
    '''
        Подсчет полной энергии системы
    '''
    if type:
        E = 0
        for part in Particles:
            E += 0.5 * MASS * (part.vx ** 2 + part.vy ** 2 + part.vz ** 2)
            for obj in Particles:
                if obj != part:
                    x, y, z = part.nearest_reflection(obj)
                    r = Distance(part, x, y, z)
                    if r != 0:
                        E += 2 * EPSILON * ((r) ** -12 - (r) ** -6)
        return E / Number_of_particles
    else:
        E = []
        for part in Particles:
            E_1 = 0.5 * MASS * (part.vx ** 2 + part.vy ** 2 + part.vz ** 2)
            for obj in Particles:
                if obj != part:
                    x, y, z = part.nearest_reflection(obj)
                    r = Distance(part, x, y, z)
                    if r != 0:
                        E_1 += 2 * EPSILON * ((r) ** -12 - (r) ** -6)
            E.append(E_1)
        return E

def K_energy(Particles):
    E = 0
    for  part in Particles:
        E += 0.5 * MASS * (part.vx ** 2 + part.vy ** 2 + part.vz ** 2)
    return E / Number_of_particles

def P_energy(Particles):
    E = 0
    for part in Particles:
        for obj in Particles:
                if obj != part:
                    x, y, z = part.nearest_reflection(obj)
                    r = Distance(part, x, y, z)
                    if r != 0:
                        E += 2 * EPSILON * ((r) ** -12 - (r) ** -6)
    return E / Number_of_particles

def Update_positions(Particles):
    '''
        Решение дифф уравнения
    '''
    Update_acceleration(Particles)
    Update_velocity(Particles)
    for part in Particles:
        x, y, z = part.x, part.y, part.z
        part.x = 2 * part.x - part.prev_x + part.ax * dt ** 2
        part.y = 2 * part.y - part.prev_y + part.ay * dt ** 2
        part.z = 2 * part.z - part.prev_z + part.az * dt ** 2
        part.prev_x, part.prev_y, part.prev_z = x, y, z
        part.x %= Size_of_box
        part.y %= Size_of_box
        part.z %= Size_of_box

def Get_velocity(Particles):
    V = []
    for part in Particles:
        V.append(np.sqrt(part.vx ** 2 + part.vy ** 2 + part.vz ** 2))
    return V

def Get_momentum(Particles):
    P = [0, 0, 0]
    for part in Particles:
        P[0] += part.vx
        P[1] += part.vy
        P[2] += part.vz
    return P

if __name__ == "__main__":
    print('This is not the main file.')