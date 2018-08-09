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

#ifndef _PARTICLE_HPP
#define _PARTICLE_HPP

#include "types.hpp"

#include <cmath>
#include <malloc.h>
#include <new>
#include <iostream>

struct Particle {

public:
    real_t *pos_x, *pos_y, *pos_z;    //position
    real_t *vel_x, *vel_y, *vel_z;    //velocity
    real_t *acc_x, *acc_y, *acc_z;    //acceleration
    real_t *mass;                     //mass

    Particle() {
        init_zero();
    }

    ~Particle() {
        dealloc();
    }

    void alloc(int nparts) {
        dealloc();

        pos_x = (real_t*) _mm_malloc(nparts * sizeof(real_t), 64);
        pos_y = (real_t*) _mm_malloc(nparts * sizeof(real_t), 64);
        pos_z = (real_t*) _mm_malloc(nparts * sizeof(real_t), 64);

        vel_x = (real_t*) _mm_malloc(nparts * sizeof(real_t), 64);
        vel_y = (real_t*) _mm_malloc(nparts * sizeof(real_t), 64);
        vel_z = (real_t*) _mm_malloc(nparts * sizeof(real_t), 64);

        acc_x = (real_t*) _mm_malloc(nparts * sizeof(real_t), 64);
        acc_y = (real_t*) _mm_malloc(nparts * sizeof(real_t), 64);
        acc_z = (real_t*) _mm_malloc(nparts * sizeof(real_t), 64);

        mass  = (real_t*) _mm_malloc(nparts * sizeof(real_t), 64);

        if ( !(pos_x && pos_y && pos_z &&
               vel_x && vel_y && vel_z &&
               acc_x && acc_y && acc_z && mass) ) {
            std::cerr << "Allocation failed (_mm_malloc returned NULL)." << std::endl;
            throw std::bad_alloc();
        }
    }

    void dealloc() {
        _mm_free(pos_x);
        _mm_free(pos_y);
        _mm_free(pos_z);

        _mm_free(vel_x);
        _mm_free(vel_y);
        _mm_free(vel_z);

        _mm_free(acc_x);
        _mm_free(acc_y);
        _mm_free(acc_z);

        _mm_free(mass);

        init_zero();
    }

private:
    void init_zero() {
        pos_x = nullptr; pos_y = nullptr; pos_z = nullptr;
        vel_x = nullptr; vel_y = nullptr; vel_z = nullptr;
        acc_x = nullptr; acc_y = nullptr; acc_z = nullptr;
        mass  = nullptr;
    }
};

#endif
