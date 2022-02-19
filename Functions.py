from ast import Num, While
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
    output = open('initial_conditions.txt', 'w')
    Particles = []
    for i in range(1, Particles_in_one_row + 1):
        for j in range(1, Particles_in_one_row + 1):
            for k in range(1, Particles_in_one_row + 1):
                Particles.append(particle(
                (i + (rn.random() - 0.5) * 2 * Rad_of_rand_gen) * Size_of_cell, 
                (j + (rn.random() - 0.5) * 2 * Rad_of_rand_gen) * Size_of_cell, 
                (k + (rn.random() - 0.5) * 2 * Rad_of_rand_gen) * Size_of_cell,
                (rn.random() - 0.5) * 2 * V_MAX,
                (rn.random() - 0.5) * 2 * V_MAX,
                (rn.random() - 0.5) * 2 * V_MAX))
                if (len(Particles) == Number_of_particles):
                    break
    for part in Particles:
        print(str(part.x) + ',' +  str(part.y) + ',' +  str(part.z) + ',' +  
        str(part.vx) + ',' +  str(part.vy) + ',' +  str(part.vz),  file=output)
    output.close()
        

def Update_acceleration(Particles):
    '''
        Производит покоординатный расчет ускорений для каждой частицы
    '''
    for part in Particles:
        part.ax = part.ay = part.az = 0
        for obj in Particles:
            if part != obj:
                x, y, z = part.nearest_reflection(obj)
                delt_x = part.x - x
                delt_y = part.y - y
                delt_z = part.z - z
                l = np.sqrt(delt_x ** 2 + delt_y ** 2 + delt_z ** 2)
                if delt_x != 0:
                    part.ax += 24 * (2 / (l ** 13) - 1 / (l ** 7)) * (delt_x / l)
                if delt_y != 0:
                    part.ay += 24 * (2 / (l ** 13) - 1 / (l ** 7)) * (delt_y / l)
                if delt_z != 0:
                    part.az += 24 * (2 / (l ** 13) - 1 / (l ** 7)) * (delt_z / l)
    
def Update_velocity(Particles):
    '''
        Рассчет моментальной скорости каждой частицы
    '''
    for part in Particles:
        part.vx += part.ax * dt
        part.vy += part.ay * dt
        part.vz += part.az * dt

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
        return E
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


def Scale(Particles):
    '''
        Попытка намутить ЗСЭ при телепортации
        Уменьшаю скорости частиц в пропорциях проекций
    '''
    E_1 = Get_energy(Particles, False)
    for part in Particles:
        part.x %= Size_of_box
        part.y %= Size_of_box
        part.z %= Size_of_box
    E_2 = Get_energy(Particles, False)
    i = 0
    for part in Particles:
        v_0 = np.sqrt((part.vx ** 2 + part.vy ** 2 + part.vz ** 2))
        dv = (E_2[i] - E_1[i]) / (MASS * v_0)
        part.vx *= 1 - dv / v_0
        part.vy *= 1 - dv / v_0
        part.vz *= 1 - dv / v_0
        i += 1

def Update_positions(Particles):
    '''
        Решение дифф уравнения методом Верле-Стёрмера
    '''
    Update_velocity(Particles)
    Update_acceleration(Particles)
    for part in Particles:
        part.x += part.vx * dt + part.ax * dt * 2
        part.y += part.vy * dt + part.ay * dt * 2
        part.z += part.vz * dt + part.az * dt * 2
        part.x %= Size_of_box
        part.y %= Size_of_box
        part.z %= Size_of_box
        #if (part.x > Size_of_box or part.x < 0 or  
        #part.y > Size_of_box or part.y < 0 or part.z > Size_of_box or part.z < 0):
            #Scale(Particles)

def Get_velocity(Particles):
    V = []
    for part in Particles:
        V.append(np.sqrt(part.vx ** 2 + part.vy ** 2 + part.vz ** 2))
    return V

def Dispersion(A):
    '''
        Проверка нормального распределения
    '''
    mean = 0
    sigma = 0
    for v in A:
        mean += np.abs(v)
    mean /= len(A)
    for v in A:
        sigma += (mean - np.abs(v)) ** 2
    sigma = np.sqrt(sigma / (len(A) - 1))
    s_1 = 0
    s_2 = 0
    s_3 = 0
    for v in A:
        if v >= mean - sigma and v <= mean + sigma:
            s_1 += 1
        if v >= mean - 2 * sigma and v <= mean + 2 * sigma:
            s_2 += 1
        if v >= mean - 3 * sigma and v <= mean + 3 * sigma:
            s_3 += 1
    print('Mean value: ' + str(mean))
    print('Dispersion: ' + str(sigma))
    print('1 sigma: ' + str(s_1 / len(A)))
    print('2 sigma: ' + str(s_2 / len(A)))
    print('3 sigma: ' + str(s_3 / len(A)))

if __name__ == "__main__":
    print('This is not the main file.')