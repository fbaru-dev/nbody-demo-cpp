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

#ifndef _GSIMULATION_HPP
#define _GSIMULATION_HPP

#include "Particle.hpp"

class GSimulation {

public:
    //Constructor and destructor
    GSimulation();
    ~GSimulation();

    //A function for simulation
    void start();

    //Public setters
    inline void set_nparts(int nparts) { _set_nparts(nparts); }
    inline void set_nsteps(int nsteps) { _set_nsteps(nsteps); }

    //Public getters
    inline int    get_nparts()   const { return _nparts;   }
    inline int    get_nsteps()   const { return _nsteps;   }
    inline real_t get_tstep()    const { return _tstep;    }
    inline real_t get_nthreads() const { return _nthreads; }

private:
    Particle* _particles;     //array of particles

    //Simulation parameters
    int       _nparts;        //number of particles
    int       _nsteps;        //number of integration steps
    real_t    _tstep;         //time step of the simulation
    int       _nthreads;      //number of threads

    //Simulation result variables
    real_t    _kenergy;       //kinetic energy

    double    _tot_time;      //total time of the simulation
    double    _tot_flops;     //total number of flops

    double    _gflops_avg;    //gflops average
    double    _gflops_dev;    //gflops deviation

    //Functions for memory operations
    void _alloc();
    void _dealloc();

    //Functions for zero initialization
    //UNSAFE! Should only be invoked by
    //constructor(), _allocate(), _deallocate()
    void _init_zero();

    //Functions for initializing particles fields with random values
    void _init_particles_pos(unsigned seed);
    void _init_particles_vel(unsigned seed);
    void _init_particles_acc();
    void _init_particles_mass(unsigned seed);

    //Private setters
    inline void _set_nparts(int nparts)  { _nparts  = nparts; }
    inline void _set_nsteps(int nsteps)  { _nsteps  = nsteps; }
    inline void _set_tstep(real_t tstep) { _tstep   = tstep;  }

    //Functions for printing report header, result lines, and footer
    void _print_header();
    void _print_resultline(int step, real_t kenergy, double time, double ngflop);
    void _print_footer();
};

#endif
