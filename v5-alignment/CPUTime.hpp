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

#ifndef _CPUTIME_HPP
#define _CPUTIME_HPP

#include <omp.h>
#include <chrono>

namespace CPUTime {

    //Implementation using OpenMP
    inline double omp_get_time_in_seconds() {
        return omp_get_wtime();
    }

    //Implementation using Chrono STD library
    inline double chrono_get_time_in_seconds() {
        namespace ch = std::chrono;
        return ch::time_point_cast<ch::milliseconds>(ch::steady_clock::now()).time_since_epoch().count() * 1.e-3;
    }

    //Implementation being used by default
    constexpr auto get_time_in_seconds = omp_get_time_in_seconds;
}

#endif
