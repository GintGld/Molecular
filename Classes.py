from math import dist
import numpy as np
from Constants import *

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
        self.prev_x = x_0 - self.vx * dt
        self.prev_y = y_0 - self.vy * dt
        self.prev_z = z_0 - self.vz * dt
import numpy as np
from Constants import *

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

if __name__ == "__main__":
    print('This is not the main file.')