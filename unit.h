#ifndef __UNIT_H__
#define __UNIT_H__

#include<iostream>
#include<iomanip>
#include<fstream>
#include<math.h>
#include<time.h>
#include<random>
#include<string>
#include<vector>
#define ld long double

using namespace std;

int number_of_particles = 3, step_of_count = 1000;
ld SIGMA = 1, EPSILON = 1, MASS = 1, time_of_relaxation = 70, dt = 0.001, temperature = 1, time_of_recording = 200,
rad_of_rand_gen = 0.1, velocity_dispersion = 1.5, concetration = 0.15, mean_temp = 0;

ld ACCELERATION, VELOCITY, TIME, size_of_cell, size_of_box, P_en, K_en;
int particles_in_one_row;

ofstream vmd, eout, vout, diff;

class particle
{
    public:
        ld x, y, z, vx, vy, vz, ax = 0, ay = 0, az = 0, prev_x, prev_y, prev_z, prev_ax, prev_ay, prev_az;
        int step_x = 0, step_y = 0, step_z = 0;
        void init_coordinates(ld x0, ld y0, ld z0)
        {
            x = x0; y = y0; z = z0;
            prev_x = x - vx * dt;
            prev_y = y - vy * dt;
            prev_z = z - vz * dt;
            return;
        }
        void init_velocity(ld vx_0, ld vy_0, ld vz_0)
        {
            vx = vx_0; vy = vy_0; vz = vz_0;
            return;
        }
        void update_acceleration(ld d_ax, ld d_ay, ld d_az)
        {
            ax += d_ax; ay += d_ay; az += d_az;
            return;
        }
};

void correct_CM(vector <particle>& particles)
{
    ld V_x = 0, V_y = 0, V_z = 0;
    for (int i = 0; i < number_of_particles; ++i) {V_x += particles[i].vx; V_y += particles[i].vy; V_z += particles[i].vz;}
    V_x /= number_of_particles; V_y /= number_of_particles; V_z /= number_of_particles;
    for (int i = 0; i < number_of_particles; ++i) {particles[i].vx -= V_x; particles[i].vy -= V_y; particles[i].vz -= V_z;}
    return;
}

void scale_velocity(vector <particle>& particles)
{
    K_en = 0;
    for (int i = 0; i < number_of_particles; ++i) {K_en += particles[i].vx * particles[i].vx + particles[i].vy * particles[i].vy + particles[i].vz * particles[i].vz;}
    K_en /= number_of_particles;
    ld k = sqrt(temperature * number_of_particles * 3/ K_en);
    for (int i = 0; i < number_of_particles; ++i) {particles[i].vx *= k; particles[i].vy *= k; particles[i].vz *= k;}
    return;
}

void set_parameters()
{
    ACCELERATION = EPSILON / (MASS * SIGMA), VELOCITY = sqrt(EPSILON / MASS), TIME = sqrt(MASS / EPSILON) * SIGMA;
    ld d = cbrt(number_of_particles) - 1. * floor(cbrt(number_of_particles));
    int f = 0;
    if (d > 0) f = 1;
    particles_in_one_row = floor(cbrt(number_of_particles)) + f;
    size_of_cell = cbrt(number_of_particles / concetration) / (particles_in_one_row + 1),
    size_of_box = (particles_in_one_row + 1) * size_of_cell;
    return;
}

ld dist(particle a, particle b)
{
    return sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y) + (a.z - b.z) * (a.z - b.z));
}

ld dist_0(particle a, ld x, ld y, ld z)
{
    return sqrt((a.x - x) * (a.x - x) + (a.y - y) * (a.y - y) + (a.z - z) * (a.z - z));
}

ld exp(ld a, int k)
{
    ld f = 1;
    for (int i = 0; i < k; ++i) f /= a;
    return f;
}

void nearest_reflection(particle a, particle b, ld c[3])
{
    ld d = dist(a, b), minim = 10 * size_of_box;
    c[0] = b.x, c[1] = b.y, c[2] = b.z;
    if (d <= 0.5 * size_of_box) return;
    for (int i = -1; i < 2; ++i) for (int j = -1; j < 2; ++j) for (int k = -1; k < 2; ++k)
    {
        d = dist_0(a, b.x + i * size_of_box, b.y + j * size_of_box, b.z + k * size_of_box);
        if (d < minim)
        {
            minim= d;
            c[0] = b.x + i * size_of_box; c[1] = b.y + j * size_of_box; c[2] = b.z + k * size_of_box;
        }
    } 
    return;
}

void generate_particles()
{
    default_random_engine generator(time(NULL));
    normal_distribution<ld> distribution(0, velocity_dispersion);
    int counter = 0;
    ofstream fout;
    fout.open("initial_conditions.txt");
    for (int i = 1; i <= particles_in_one_row; ++i) for (int j = 1; j <= particles_in_one_row; ++j) for (int k = 1; k <= particles_in_one_row; ++k)
    {
        if (counter == number_of_particles) 
        {
            cout << "Particles were generated." << endl;
            return;
        }
        fout << (i + (2 * (ld)(rand()) / RAND_MAX - 1) * rad_of_rand_gen) * size_of_cell << ',' << (j + (2 * (ld)(rand()) / RAND_MAX - 1) * rad_of_rand_gen) * size_of_cell << ',' << 
        (k + (2 * (ld)(rand()) / RAND_MAX - 1) * rad_of_rand_gen) * size_of_cell << ',' << 
        distribution(generator) << ',' << distribution(generator) << ',' << distribution(generator) << endl;
        counter++;
    }
    fout.close();
    cout << "Particles were generated." << endl;
    return;
}

void read_particles(vector <particle>& particles)
{
    ifstream fin;
    fin.open("initial_conditions.txt");
    string s;
    ld d[6];
    particle p;
    for (int i = 0; i < number_of_particles; ++i)
    {
        getline(fin, s);
        for (int i = 0; i < 6; ++i)
        {
            d[i] = stold(s.substr(0, s.find(',')));
            if (i < 5) s = s.substr(s.find(',') + 1, s.size() - s.find(','));
        }
        p.init_velocity(d[3], d[4], d[5]);
        p.init_coordinates(d[0], d[1], d[2]);
        particles.push_back(p);
    }
    fin.close();
    return;
}

void reset_files()
{
    vmd.open("Coordinates_VMD.txt"); vmd.close();
    eout.open("Energy.txt"); 
    eout << dt << "\n";
    eout.close();
    vout.open("velocity.txt"); vout.close();
    diff.open("diffusion.txt");
    for (int i = 0; i < number_of_particles; ++i) diff << i << "x," << i << "y," << i << "z,";
    diff << "N,dt\n";
    //diff << number_of_particles << "," << dt << "," << time_of_recording / dt << "\n";
    diff.close();
    return;
}

void close_files(vector <particle>& particles)
{
    eout.open("Energy.txt", ios::app); eout << "end\n"; eout.close();
    vout.open("velocity.txt", ios::app); vout << "end\n"; vout.close();
    diff.open("diffusion.txt", ios::app); diff << "end\n"; diff.close();
    return;
}

void update_acceleration(vector <particle>& particles)
{
    ld a_x, a_y, a_z, L, k;
    P_en = 0;
    ld *c = new ld [3];
    for (int i = 0; i < number_of_particles; ++i)
    {
        particles[i].prev_ax = particles[i].ax; particles[i].ax = 0;
        particles[i].prev_ay = particles[i].ay; particles[i].ay = 0;
        particles[i].prev_az = particles[i].az; particles[i].az = 0;
    }
    for (int i = 0; i < number_of_particles; ++i)
    {
        for (int j = i + 1; j < number_of_particles; ++j)
        {
            nearest_reflection(particles[i], particles[j], c);
            L = dist_0(particles[i], c[0], c[1], c[2]);
            if (true)//(L < 2.5)
            {
                k = 24 * (2 * exp(L, 14) - exp(L, 8));
                a_x = (particles[i].x - c[0]) * k;
                a_y = (particles[i].y - c[1]) * k;
                a_z = (particles[i].z - c[2]) * k;
                particles[i].update_acceleration(a_x, a_y, a_z);
                particles[j].update_acceleration(-a_x, -a_y, -a_z);
            }
            P_en += 4 * (exp(L, 12) - exp(L, 6));
        }
        if (particles[i].prev_ax == 0 && particles[i].prev_ay == 0 && particles[i].prev_az == 0)
        {
            particles[i].prev_ax = particles[i].ax; particles[i].prev_ay = particles[i].ay; particles[i].prev_az = particles[i].az;
        }
    }
    delete[] c;
    return;
}

void update_position(vector <particle>& particles)
{
    update_acceleration(particles);
    ld x, y, z;
    for (int i = 0; i < number_of_particles; ++i)
    {
        x = particles[i].x; y = particles[i].y; z = particles[i].z;
        particles[i].x += particles[i].x - particles[i].prev_x + particles[i].ax * dt * dt;
        particles[i].y += particles[i].y - particles[i].prev_y + particles[i].ay * dt * dt;
        particles[i].z += particles[i].z - particles[i].prev_z + particles[i].az * dt * dt;
        particles[i].prev_x = x; particles[i].prev_y = y; particles[i].prev_z = z;
        particles[i].vx += 0.5 * dt * (particles[i].prev_ax + particles[i].ax);
        particles[i].vy += 0.5 * dt * (particles[i].prev_ay + particles[i].ay);
        particles[i].vz += 0.5 * dt * (particles[i].prev_az + particles[i].az);
        if (particles[i].x > size_of_box) {particles[i].x -= size_of_box; ++particles[i].step_x;} // при рассчете диффузии определяем зависимость среднеквадратического смещения от температуры
        if (particles[i].y > size_of_box) {particles[i].y -= size_of_box; ++particles[i].step_y;} // (сделано) при переходе через границу надо сделать счетчик, что считать смещение по-честному
        if (particles[i].z > size_of_box) {particles[i].z -= size_of_box; ++particles[i].step_z;}
        if (particles[i].x < 0) {particles[i].x += size_of_box; --particles[i].step_x;}
        if (particles[i].y < 0) {particles[i].y += size_of_box; --particles[i].step_y;}
        if (particles[i].z < 0) {particles[i].z += size_of_box; --particles[i].step_z;}
    }
    K_en = 0, P_en = 0;
    for (int i = 0; i < number_of_particles; ++i)
    {
        ld d;
        ld *c = new ld [3];
        for(int j = i + 1; j < number_of_particles; ++j)
        {
            nearest_reflection(particles[i], particles[j], c);
            d = dist_0(particles[i], c[0], c[1], c[2]);
            P_en += 4 * (exp(d, 12) - exp(d, 6));
        }
        K_en += 0.5 * (particles[i].vx * particles[i].vx + particles[i].vy * particles[i].vy + particles[i].vz * particles[i].vz);
        delete[] c;
    }
    K_en /= number_of_particles; P_en /= number_of_particles;
    return;
}

void vmd_out(vector <particle>& particles)
{
    vmd.open("Coordinates_VMD.txt", ios::app);
    vmd << setprecision(7) << number_of_particles << endl << endl;
    for (int i = 0; i < number_of_particles; ++i) vmd << setprecision(7) << "1 " << particles[i].x << ' ' << particles[i].y << ' ' << particles[i].z << endl;
    vmd.close();
    return;
}

void energy()
{
    eout.open("Energy.txt", ios::app);
    eout << setprecision(7) << P_en << ',' << K_en << endl;
    eout.close();
    return;
}

void velocity_distribution(vector <particle>& particles)
{
    vout.open("velocity.txt", ios::app);
    for (int i = 0; i < number_of_particles; ++i) vout << particles[i].vx << ',' << particles[i].vy << ',' << particles[i].vz << endl;
    vout.close();
    return;
}

void diffusion(vector <particle>& particles)
{
    // Усреднить квадрат смещения. Зависимость от времени
    diff.open("diffusion.txt", ios::app);
    for (int i = 0; i < number_of_particles - 1; ++i)
    {
        diff << particles[i].x + particles[i].step_x * size_of_box << ',' <<
                particles[i].y + particles[i].step_y * size_of_box << ',' <<
                particles[i].z + particles[i].step_z * size_of_box << ',';
    }
    diff << number_of_particles << ',' << dt << "\n";
    diff.close();
}

void do_modeling(vector <particle>& particles)
{
    for (int step = 1; (step - 1) * dt <= time_of_relaxation + time_of_recording; ++step)
    {
        vmd_out(particles);
        update_position(particles);
        energy();
        if (step % step_of_count == 0) cout << step * dt << " time were modeled...\n";
        if (step * dt > time_of_relaxation) {velocity_distribution(particles); diffusion(particles);}
        mean_temp = K_en;
    }
}

#endif