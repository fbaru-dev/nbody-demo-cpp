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
#include <sdlt/sdlt.h>

struct Particle {

public:
    real_t pos_x, pos_y, pos_z;    //position
    real_t vel_x, vel_y, vel_z;    //velocity
    real_t acc_x, acc_y, acc_z;    //acceleration
    real_t mass;                   //mass

    Particle() {
        init_zero();
    }

    void init_zero() {
        pos_x = 0.; pos_y = 0.; pos_z = 0.;
        vel_x = 0.; vel_y = 0.; vel_z = 0.;
        acc_x = 0.; acc_y = 0.; acc_z = 0.;
        mass  = 0.;
    }
};

SDLT_PRIMITIVE(Particle, pos_x, pos_y, pos_z, vel_x, vel_y, vel_z, acc_x, acc_y, acc_z, mass)

#endif
