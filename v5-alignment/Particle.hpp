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

        pos_x = new real_t[nparts];
        pos_y = new real_t[nparts];
        pos_z = new real_t[nparts];

        vel_x = new real_t[nparts];
        vel_y = new real_t[nparts];
        vel_z = new real_t[nparts];

        acc_x = new real_t[nparts];
        acc_y = new real_t[nparts];
        acc_z = new real_t[nparts];

        mass  = new real_t[nparts];
    }

    void dealloc() {
        delete[] pos_x;
        delete[] pos_y;
        delete[] pos_z;

        delete[] vel_x;
        delete[] vel_y;
        delete[] vel_z;

        delete[] acc_x;
        delete[] acc_y;
        delete[] acc_z;

        delete[] mass;

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
