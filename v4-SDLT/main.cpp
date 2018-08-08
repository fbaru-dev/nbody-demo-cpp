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

#include <string>
#include <iostream>

int main(int argc, char** argv) {
    GSimulation sim;

    try {
        //Checks CLI arguments
        if (argc > 1) {

            //Sets number of particles if the first argument is castable to int
            sim.set_nparts(std::stoi(argv[1]));

            if (argc > 2) {
                //Sets number of steps if the second argument is castable to int
                sim.set_nsteps(std::stoi(argv[2]));
            }
        }
    } catch (const std::invalid_argument& e) {
        std::cerr << "Invalid argument. Default value will be used." << std::endl;
    } catch (const std::out_of_range& e) {
        std::cerr << "Out of range argument. Default value will be used." << std::endl;
    }

    sim.start();

    return 0;
}
