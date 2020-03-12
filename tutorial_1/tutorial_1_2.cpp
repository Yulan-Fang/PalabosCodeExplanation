/* This file is part of the Palabos library.
 *
 * The Palabos softare is developed since 2011 by FlowKit-Numeca Group Sarl
 * (Switzerland) and the University of Geneva (Switzerland), which jointly
 * own the IP rights for most of the code base. Since October 2019, the
 * Palabos project is maintained by the University of Geneva and accepts
 * source code contributions from the community.
 * 
 * Contact:
 * Jonas Latt
 * Computer Science Department
 * University of Geneva
 * 7 Route de Drize
 * 1227 Carouge, Switzerland
 * jonas.latt@unige.ch
 *
 * The most recent release of Palabos can be downloaded at 
 * <https://palabos.unige.ch/>
 *
 * The library Palabos is free software: you can redistribute it and/or
 * modify it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * The library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/* Code 1.2 in the Palabos tutorial
 */

#include "palabos2D.h"
#include "palabos2D.hh"
#include <iostream>
#include <iomanip>

using namespace plb;
using namespace std;

typedef double T;
#define DESCRIPTOR plb::descriptors::D2Q9Descriptor

const plint maxIter = 1000; // Iterate during 1000 steps.
const plint nx = 600;       // Choice of lattice dimensions.
const plint ny = 600;
const T omega = 1.;        // Choice of the relaxation parameter

T rho0 = 1.; // All cells have initially density rho ...
// .. except for those inside the disk which have density
//    rho+deltaRho
T deltaRho = 1.e-4;
Array<T,2> u0(0,0);
/*Chinese version:
*基本的代码解释都在tutorial_1_1，如非必要不重复解释。
*tutorial_1_1里，你会发现"void xxxx(MultiBlockLattice2D<T,DESCRIPTOR>& lattice)"。
*在这个功能里面首先通过用lattice.getNx()和lattice.getNy()得到nx,ny来组成函数让Box2D来定义域的范围。
*随后用initializeAtEquilibrium在定义rho和u。
*末尾是lattice.innitialize()，并且在程序后面的主循环处还会有一行“xxxx(lattice);”。
*在本算例下面有一个 "void xxxx(plint iX, plint iY, T& rho, Array<T,2>& u)"。
*82-85行与tutorial_1_1里面设定常数rho的效果一样
*通过这种方法我们可以先定义子域范围和rho的值，然后在"initializeAtEquilibrium"里使用，如第110行。
*/
/* English version:
* The basic explanations are in tutorial_1_1. There will be only a little repeats of the explanation.
* In tutorial_1_1 we see function like "void xxxx(MultiBlockLattice2D<T,DESCRIPTOR>& lattice)".
* In which the area was firstly defined by Box2D by function of nx, ny from lattice.getNx(), lattice.getNy().
* Then initializeAtEquilibrium was used for difining the rho and u.
* With the end of lattice.innitialize(). Also before the main loop should we add a line of "xxxx(lattice);".
* Here we see "void xxxx(plint iX, plint iY, T& rho, Array<T,2>& u)".
* Line 82-85 has the same effect as in the tutorial_1_1 for defining the the constant rho in a sub-domain.
* and by this way we define rho to sub-domains for use inside "initializeAtEquilibrium", for example in line 110.
*
* Code explanation by Yulan Fang 
* Error correction please send to ahdhfang@hotmail.com
* ——March,12 2020 at Siwa, Egypt.
*/
void initializeConstRho(plint iX, plint iY, T& rho, Array<T,2>& u) {
    u = u0;
    rho = rho0 + deltaRho;
}

void initializeRhoOnDisk(plint iX, plint iY, T& rho, Array<T,2>& u) {
    plint radius = nx/6;
    plint centerX = nx/3;
    plint centerY = ny/4;
    u = u0;
    if( (iX-centerX)*(iX-centerX) + (iY-centerY)*(iY-centerY) < radius*radius) {
        rho = rho0 + deltaRho;
    }
    else {
        rho = rho0;
    }
}

// Initialize the lattice at zero velocity and constant density, except
//   for a slight density excess on a circular sub-domain.
void defineInitialDensityAtCenter(MultiBlockLattice2D<T,DESCRIPTOR>& lattice)
{
    // Initialize constant density everywhere.
    initializeAtEquilibrium (
           lattice, lattice.getBoundingBox(), rho0, u0 );

    // And slightly higher density in the central box.
    initializeAtEquilibrium (
           lattice, lattice.getBoundingBox(), initializeRhoOnDisk );

    lattice.initialize();
}

int main(int argc, char* argv[]) {
    plbInit(&argc, &argv);
    global::directories().setOutputDir("./tmp/");

    MultiBlockLattice2D<T, DESCRIPTOR> lattice (
           nx, ny, new BGKdynamics<T,DESCRIPTOR>(omega) );

    lattice.periodicity().toggleAll(true); // Set periodic boundaries.

    defineInitialDensityAtCenter(lattice);

    // Main loop over time iterations.
    for (plint iT=0; iT<maxIter; ++iT) {
        if (iT%40==0) {  // Write an image every 40th time step.
            pcout << "Writing GIF file at iT=" << iT << endl;
            // Instantiate an image writer with the color map "leeloo".
            ImageWriter<T> imageWriter("leeloo");
            // Write a GIF file with colors rescaled to the range of values
            //   in the matrix
            imageWriter.writeScaledGif (
                    createFileName("u", iT, 6),
                    *computeVelocityNorm(lattice) );
        }
        // Execute lattice Boltzmann iteration.
        lattice.collideAndStream();
    }
}
