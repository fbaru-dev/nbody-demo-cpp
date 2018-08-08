/****************************************************************************

    This file is part of the example codes which have been used
    for the "Code Optimization Workshop".

    Copyright (C) 2016, 2018

    Original Author: Fabio Baruffa <fbaru-dev@gmail.com>
    Modified by:  Aleksei Abliazov <aleksei.abliazov@gmail.com>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.

****************************************************************************/

#include "GSimulation.hpp"
#include "CPUTime.hpp"

#include <random>
#include <iomanip>
#include <iostream>

//----------------------------//
// Constructor and destructor //
//----------------------------//

GSimulation::GSimulation() {
    _set_nparts(16000);
    _set_nsteps(10);
    _set_tstep(0.1);

     //Set the pointers to null
     //NB:
     //we will allocate/deallocate memory just before/after the simulation
     //always keeping the pointers null when the memory is not allocated
     //(initialize them with null after every deallocation and initially in constructor)
    _init_zero();
}

GSimulation::~GSimulation() {
    _dealloc();                                 //deallocate memory and set ptrs to null
}

//---------------------------------//
// Functions for memory operations //
//---------------------------------//

void GSimulation::_alloc() {
    _dealloc();                                 //deallocate memory and set ptrs to null
    _particles = new Particle();                //allocate memory
    _particles->alloc(get_nparts());
}

void GSimulation::_dealloc() {
    delete _particles;                          //deallocate memory
    _init_zero();                               //set pointers to null
}

//-------------------------------------------//
// Functions for zero initialization         //
// UNSAFE! Should only be invoked by         //
// constructor(), _allocate(), _deallocate() //
//-------------------------------------------//

void GSimulation::_init_zero() {
    _particles = nullptr;                       //set pointers to null
}

//----------------------------------------------------------------//
// Functions for initializing particles fields with random values //
//----------------------------------------------------------------//

void GSimulation::_init_particles_pos(unsigned seed) {
    int nparts = get_nparts();

    std::mt19937 gen(seed);
    std::uniform_real_distribution<real_t> unif_d(0, 1.0);

    for (int i = 0; i < nparts; ++i) {
        _particles->pos_x[i] = unif_d(gen);
        _particles->pos_y[i] = unif_d(gen);
        _particles->pos_z[i] = unif_d(gen);
    }
}

void GSimulation::_init_particles_vel(unsigned seed) {
    int nparts = get_nparts();

    std::mt19937 gen(seed);
    std::uniform_real_distribution<real_t> unif_d(-1.0, 1.0);

    for (int i = 0; i < nparts; ++i) {
        _particles->vel_x[i] = unif_d(gen) * 1.0e-3f;
        _particles->vel_y[i] = unif_d(gen) * 1.0e-3f;
        _particles->vel_z[i] = unif_d(gen) * 1.0e-3f;
    }
}

void GSimulation::_init_particles_acc() {
    int nparts = get_nparts();

    for (int i = 0; i < nparts; ++i) {
        _particles->acc_x[i] = 0.f;
        _particles->acc_y[i] = 0.f;
        _particles->acc_z[i] = 0.f;
    }
}

void GSimulation::_init_particles_mass(unsigned seed) {
    int nparts = get_nparts();

    std::mt19937 gen(seed);
    std::uniform_real_distribution<real_t> unif_d(0.0, 1.0);

    real_t n = static_cast<real_t>(nparts);

    for (int i = 0; i < nparts; ++i) {
        _particles->mass[i] = unif_d(gen) * n;
    }
}

//----------------------------------------------------------------//
// Functions for printing report header, result lines, and footer //
//----------------------------------------------------------------//

void GSimulation::_print_header() {
    std::cout << "================================================" << std::endl
              << "# Gravity Simulation"                             << std::endl;

    std::cout << "# nPart: "  << get_nparts()
              << ", nSteps: " << get_nsteps()
              << ", dt: "     << get_tstep()
              << std::endl    << std::endl;

    std::cout << "------------------------------------------------" << std::endl;

    std::cout << " "
              << std::left << std::setw(8)  << "s"
              << std::left << std::setw(8)  << "dt"
              << std::left << std::setw(12) << "kenergy"
              << std::left << std::setw(12) << "time (s)"
              << std::left << std::setw(12) << "GFlops"
              << std::endl;

    std::cout << "------------------------------------------------" << std::endl;
}

void GSimulation::_print_resultline(int step, real_t kenergy, double time, double ngflop) {
    std::cout << " "
              << std::left << std::setw(8)                          << step
              << std::left << std::setprecision(5) << std::setw(8)  << step * get_tstep()
              << std::left << std::setprecision(5) << std::setw(12) << kenergy
              << std::left << std::setprecision(5) << std::setw(12) << time
              << std::left << std::setprecision(5) << std::setw(12) << ngflop / time
              << std::endl;
}

void GSimulation::_print_footer() {
    std::cout << "------------------------------------------------" << std::endl;

    std::cout                                                       << std::endl
              << "# Number Threads   : " << _nthreads               << std::endl
              << "# Total Time (s)   : " << _tot_time               << std::endl
              << "# Avg. Performance : " << _gflops_avg
                                         << " +- "
                                         << _gflops_dev             << std::endl;

    std::cout << "================================================" << std::endl;
}

//---------------------------//
// A function for simulation //
//---------------------------//

void GSimulation::start() {

    //------------------------------------------//
    // Variables declaration and initialization //
    //------------------------------------------//

    //Simulation independent constants
    const real_t G = 6.67259e-11;            //gravitational constant
    const real_t softeningSquared = 1e-3;    //prevents explosion if the particles
                                             //are close to each other

    //Simulation dependent constants
    const int    nparts = get_nparts();      //number of particles
    const int    nsteps = get_nsteps();      //number of integration steps
    const real_t dt     = get_tstep();       //time step of the simulation


    //Simulation result variables
    real_t step_kenergy        = 0.;         //kinetic energy in each step
    double step_time_start     = 0.;         //step starting time
    double step_time_duration  = 0.;         //step duration
    double total_time_start    = 0.;         //simulation starting time
    double total_time_duration = 0.;         //simulation duration

    double gflops_avg          = 0.;         //gflops average
    double gflops_dev          = 0.;         //glofps deviation


    //Temporary variables
    int    step, i, j;                       //vars used for iteration
    real_t dx, dy, dz;                       //xyz distance
    real_t distanceSqr;                      //squared distance
    real_t distanceInv;                      //1/distance


    //Number of gflop in each step
    //  1e-9          : flop -> gflop
    //  (11. + 18.)   : number of flop for each iteration (the first loop)
    //  npart * npart : number of iterations              (the first loop)
    //  (19.)         : number of flop for each iteration (the second loop)
    //  npart         : number of iterations              (the second loop)
    double ngflop = 1e-9 * ((11. + 18.) * nparts * nparts + (19.) * nparts);


    //-----------------------------------------//
    // Particles allocation and initialization //
    //-----------------------------------------//

    _alloc();                                //allocates memory

    //std::random_device rd;                 //you can use
    const unsigned seed = 42;                //rd(); as a seed

    _init_particles_pos(seed);               //initializes position
    _init_particles_vel(seed);               //initializes velocity
    _init_particles_acc();                   //initializes acceleration
    _init_particles_mass(seed);              //initializes mass


    //------------------//
    // Simulation start //
    //------------------//

    _print_header();

    total_time_start = CPUTime::get_time_in_seconds();

    //In each step
    for (step = 1; step <= nsteps; ++step) {

        step_time_start = CPUTime::get_time_in_seconds();

        //Iterates over all particles
        for (i = 0; i < nparts; ++i) {

            //Resets acceleration
            _particles->acc_x[i] = 0.f;
            _particles->acc_y[i] = 0.f;
            _particles->acc_z[i] = 0.f;

            //For given particle
            //computes the distance to other particles
            //and updates acceleration
            //using Newton's law of gravitation
            for (j = 0; j < nparts; ++j) {

                //Computes the distance
                dx = _particles->pos_x[j] - _particles->pos_x[i];                   //1flop
                dy = _particles->pos_y[j] - _particles->pos_y[i];                   //1flop
                dz = _particles->pos_z[j] - _particles->pos_z[i];                   //1flop

                distanceSqr = dx * dx + dy * dy + dz * dz + softeningSquared;       //6flops
                distanceInv = 1.0f / sqrtf(distanceSqr);                            //1div+1sqrt

                //Updates acceleration
                _particles->acc_x[i] += G * _particles->mass[j] * dx *
                                        distanceInv * distanceInv * distanceInv;    //6flops
                _particles->acc_y[i] += G * _particles->mass[j] * dy *
                                        distanceInv * distanceInv * distanceInv;    //6flops
                _particles->acc_z[i] += G * _particles->mass[j] * dz *
                                        distanceInv * distanceInv * distanceInv;    //6flops
            }
        }

        //Resets kinetic energy for given iteration step
        step_kenergy = 0;

        //Iterates over all particles
        for (i = 0; i < nparts; ++i) {

            //Updates velocity for given particle
            _particles->vel_x[i] += _particles->acc_x[i] * dt;                      //2flops
            _particles->vel_y[i] += _particles->acc_y[i] * dt;                      //2flops
            _particles->vel_z[i] += _particles->acc_z[i] * dt;                      //2flops

            //Updates position for given particle
            _particles->pos_x[i] += _particles->vel_x[i] * dt;                      //2flops
            _particles->pos_y[i] += _particles->vel_y[i] * dt;                      //2flops
            _particles->pos_z[i] += _particles->vel_z[i] * dt;                      //2flops

            //Adds particle kinetic energy to step kinetic energy                   //7flops
            step_kenergy += _particles->mass[i] * (_particles->vel_x[i] * _particles->vel_x[i] +
                                                   _particles->vel_y[i] * _particles->vel_y[i] +
                                                   _particles->vel_z[i] * _particles->vel_z[i] );
        }

        //Kinetic energy at the current step
        _kenergy = 0.5 * step_kenergy;

        //Step duration in seconds for the current step
        step_time_duration = CPUTime::get_time_in_seconds() - step_time_start;

        //Computes the average and deviation
        gflops_avg += ngflop / step_time_duration;
        gflops_dev += ngflop * ngflop / (step_time_duration * step_time_duration);

        //Prints the results line
        _print_resultline(step, _kenergy, step_time_duration, ngflop);
    }

    //Total duration in seconds for the simulation
    total_time_duration = CPUTime::get_time_in_seconds() - total_time_start;

    //Total time and gflops
    _tot_time    = total_time_duration;
    _tot_flops   = ngflop * nsteps;

    _nthreads = 1;

    _gflops_avg =      gflops_avg / (double) (nsteps);
    _gflops_dev = sqrt(gflops_dev / (double) (nsteps) - _gflops_avg * _gflops_avg);

    //Prints the footer
    _print_footer();


    //------------------------//
    // Particles deallocation //
    //------------------------//

    _dealloc();
}
