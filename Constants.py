import numpy as np

SIGMA = 1
EPSILON = 1
MASS = 1
ACCELERATION = EPSILON / (MASS * SIGMA)
VELOCITY = np.sqrt(EPSILON / MASS)
TIME = np.sqrt(MASS / EPSILON) * SIGMA

Number_of_particles = 10
d = Number_of_particles ** (1 / 3) - int(Number_of_particles ** (1 / 3))
if (d > 0): d = 1
Particles_in_one_row = int(Number_of_particles ** (1 / 3)) + d
Concetration = 0.1 * SIGMA ** (-3)
Size_of_cell = Concetration ** (-1 / 3)
Size_of_box = (Particles_in_one_row + 1) * Size_of_cell
Rad_of_rand_gen = 0.2
Number_of_steps = 1500
Step_of_count = 500

V_MAX = 1.5
dt = 0.001

if __name__ == "__main__":
    #print('This is not the main file.')
    print(Size_of_cell)