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
}

GSimulation::~GSimulation() {
}

//---------------------------------//
// Functions for memory operations //
//---------------------------------//

void GSimulation::_alloc() {
    _particles.resize(get_nparts());
}

void GSimulation::_dealloc() {
    //RAII, So no need to deallocate manually
}

//-------------------------------------------//
// Functions for zero initialization         //
// UNSAFE! Should only be invoked by         //
// constructor(), _allocate(), _deallocate() //
//-------------------------------------------//

void GSimulation::_init_zero() {
}

//----------------------------------------------------------------//
// Functions for initializing particles fields with random values //
//----------------------------------------------------------------//

void GSimulation::_init_particles_pos(unsigned seed) {
    int nparts = get_nparts();

    std::mt19937 gen(seed);
    std::uniform_real_distribution<real_t> unif_d(0, 1.0);

    auto particles = _particles.access();
    for (int i = 0; i < nparts; ++i) {
        particles[i].pos_x() = unif_d(gen);
        particles[i].pos_y() = unif_d(gen);
        particles[i].pos_z() = unif_d(gen);
    }
}

void GSimulation::_init_particles_vel(unsigned seed) {
    int nparts = get_nparts();

    std::mt19937 gen(seed);
    std::uniform_real_distribution<real_t> unif_d(-1.0, 1.0);

    auto particles = _particles.access();
    for (int i = 0; i < nparts; ++i) {
        particles[i].vel_x() = unif_d(gen) * 1.0e-3f;
        particles[i].vel_y() = unif_d(gen) * 1.0e-3f;
        particles[i].vel_z() = unif_d(gen) * 1.0e-3f;
    }
}

void GSimulation::_init_particles_acc() {
    int nparts = get_nparts();

    auto particles = _particles.access();
    for (int i = 0; i < nparts; ++i) {
        particles[i].acc_x() = 0.f;
        particles[i].acc_y() = 0.f;
        particles[i].acc_z() = 0.f;
    }
}

void GSimulation::_init_particles_mass(unsigned seed) {
    int nparts = get_nparts();

    std::mt19937 gen(seed);
    std::uniform_real_distribution<real_t> unif_d(0.0, 1.0);

    real_t n = static_cast<real_t>(nparts);

    auto particles = _particles.access();
    for (int i = 0; i < nparts; ++i) {
        particles[i].mass() = unif_d(gen) * n;
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

    auto particles = _particles.access();

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
            real_t acc_x = 0.f;
            real_t acc_y = 0.f;
            real_t acc_z = 0.f;

            //For given particle
            //computes the distance to other particles
            //and updates acceleration
            //using Newton's law of gravitation
            for (j = 0; j < nparts; ++j) {

                //Computes the distance
                dx = particles[j].pos_x() - particles[i].pos_x();                   //1flop
                dy = particles[j].pos_y() - particles[i].pos_y();                   //1flop
                dz = particles[j].pos_z() - particles[i].pos_z();                   //1flop

                distanceSqr = dx * dx + dy * dy + dz * dz + softeningSquared;       //6flops
                distanceInv = 1.0f / sqrtf(distanceSqr);                            //1div+1sqrt

                //Updates acceleration
                acc_x += G * particles[j].mass() * dx *
                         distanceInv * distanceInv * distanceInv;                   //6flops
                acc_y += G * particles[j].mass() * dy *
                         distanceInv * distanceInv * distanceInv;                   //6flops
                acc_z += G * particles[j].mass() * dz *
                         distanceInv * distanceInv * distanceInv;                   //6flops
            }

            particles[i].acc_x() = acc_x;
            particles[i].acc_y() = acc_y;
            particles[i].acc_z() = acc_z;
        }

        //Resets kinetic energy for given iteration step
        step_kenergy = 0;

        //Iterates over all particles
        #pragma ivdep
        for (i = 0; i < nparts; ++i) {

            //Updates velocity for given particle
            particles[i].vel_x() += particles[i].acc_x() * dt;                      //2flops
            particles[i].vel_y() += particles[i].acc_y() * dt;                      //2flops
            particles[i].vel_z() += particles[i].acc_z() * dt;                      //2flops

            //Updates position for given particle
            particles[i].pos_x() += particles[i].vel_x() * dt;                      //2flops
            particles[i].pos_y() += particles[i].vel_y() * dt;                      //2flops
            particles[i].pos_z() += particles[i].vel_z() * dt;                      //2flops

            //Adds particle kinetic energy to step kinetic energy                   //7flops
            step_kenergy += particles[i].mass() * (particles[i].vel_x() * particles[i].vel_x() +
                                                   particles[i].vel_y() * particles[i].vel_y() +
                                                   particles[i].vel_z() * particles[i].vel_z() );
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
