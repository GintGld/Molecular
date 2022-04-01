#include<iostream>
#include<fstream>
#include<math.h>
#include<time.h>
#include<random>
#include<string>
#include<vector>
#include<iomanip>
//#include<sstream>
#define ld long double

#pragma hdrstop

#include<c:\\users\\coolg\\molecular\\unit.h>

using namespace std;

int main()
{
    srand(time(NULL));
    set_parameters();
    vector <particle> particles;
    generate_particles();
    read_particles(particles);
    correct_CM(particles);

    reset_files();
    do_modeling(particles);
    close_files(particles);
    cout << "Modeling ended.";
    cout << "\nTemp: " << mean_temp * 2 / 3;
    return 0;
}
