import matplotlib.pyplot as plt
import numpy as np
import random as rn

SIGMA = 1
EPSILON = 1
MASS = 1
Number_of_particles = 50
Concetration = 0.1

ACCELERATION = EPSILON / (MASS * SIGMA)
VELOCITY = np.sqrt(EPSILON / MASS)
TIME = np.sqrt(MASS / EPSILON) * SIGMA

d = Number_of_particles ** (1 / 3) - int(Number_of_particles ** (1 / 3))
if (d > 0): d = 1
Particles_in_one_row = int(Number_of_particles ** (1 / 3)) + d
Size_of_cell = (Number_of_particles / Concetration) ** (1 / 3) / (Particles_in_one_row + 1)
Size_of_box = (Particles_in_one_row + 1) * Size_of_cell
Rad_of_rand_gen = 0.2
Number_of_steps = 100
Step_of_count = 500
Velocity_dispersion = 0.5
dt = 0.001

Particles = []

#   Задание класса

def Distance(object, x, y, z):
    return np.sqrt((object.x - x) ** 2 + (object.y - y) ** 2 + (object.z - z) ** 2)

class particle:
    def __init__(self, x_0, y_0, z_0, v_x, v_y, v_z):
        '''
            Создание объекта частица с заданными координатами и скоростями
        '''
        self.x = x_0
        self.y = y_0
        self.z = z_0
        self.vx = v_x
        self.vy = v_y
        self.vz = v_z
        self.prev_x = self.x - self.vx * dt
        self.prev_y = self.y - self.vy * dt
        self.prev_z = self.z - self.vz * dt
        self.prev_ax = 0
        self.prev_ay = 0
        self.prev_az = 0
        self.ax = 0
        self.ay = 0
        self.az = 0

    def nearest_reflection(self, object):
        '''
            Возвращает ближайший образ точки (object) к данной (self)
        '''
        x = object.x
        y = object.y
        z = object.z
        x_0 = x
        y_0 = y
        z_0 = z
        L = Size_of_box
        minim = 10 * Size_of_box
        if Distance(self, x, y, z) <= 0.5 * Size_of_box:
            return x, y, z
        for i in [-1, 0, 1]:
            for j in [-1, 0, 1]:
                for k in [-1, 0, 1]:
                    if minim > Distance(self, x + i * L, y + j * L, z + k * L):
                        x_0 = x + i * L
                        y_0 = y + j * L
                        z_0 = z + k * L
                        minim = Distance(self, x + i * L, y + j * L, z + k * L)
        return x_0, y_0, z_0

#   Вспомогательные функции для моделирования и отображения результатов

def Generate_particles():
    '''
        Генерирует массив частиц, расположенных на расстоянии не больших четверти шага сетки 
        с рандомными  скоростями, ограниченные константой  V_MAX
        Возвращает массив с элементами класса particle
    '''
    vx = np.random.normal(0, Velocity_dispersion, Number_of_particles)
    vy = np.random.normal(0, Velocity_dispersion, Number_of_particles)
    vz = np.random.normal(0, Velocity_dispersion, Number_of_particles)
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

#   Основные функции, вызываются в main()

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

    print(str(Number_of_steps) + ' of ' + str(Number_of_steps) + ' were modeled...')
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
    Momentum()