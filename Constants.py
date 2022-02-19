import numpy as np

SIGMA = 1
EPSILON = 1
Size_of_cell = 4 * SIGMA
MASS = 1
ACCELERATION = EPSILON / (MASS * SIGMA)
VELOCITY = np.sqrt(EPSILON / MASS)
TIME = np.sqrt(MASS / EPSILON) * SIGMA

Particles_in_one_row = 3
Number_of_particles = 15
Size_of_box = (Particles_in_one_row + 1) * Size_of_cell
Rad_of_rand_gen = 0.2
Number_of_steps = 20000
Step_of_count = 1000

V_MAX = 5
dt = 0.001

if __name__ == "__main__":
    print('This is not the main file.')