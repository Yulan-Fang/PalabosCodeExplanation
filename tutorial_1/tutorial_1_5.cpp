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

#include "palabos2D.h"
#include "palabos2D.hh"
#include <vector>
#include <iostream>
#include <iomanip>

/* Code 1.5 in the Palabos tutorial
 */

using namespace plb;
using namespace std;

typedef double T;
#define DESCRIPTOR plb::descriptors::D2Q9Descriptor

/// Velocity on the parabolic Poiseuille profile
//本tutorial主要是介绍通过类来定义域和IncomprFlowParam储存计算参数
//此处计算得出抛物线形式的泊肃叶速度场线，其中中间位置速度值最高
//初见你可能会觉得一头雾水，这个IncomprFlowParam是个什么东西，不过后面有说明
/*In this tutorial file you will see how to define a domain by a class and how IncomprFlowParam stores parameters.
*Below it is calculation of Poisueille profile, and the middle of y-axis has the highest velocity value. 
*You may get puzzled by this IncomprFlowParam, later it will be explained.
*/
T poiseuilleVelocity(plint iY, IncomprFlowParam<T> const& parameters) {
    T y = (T)iY / parameters.getResolution();
    //这里他y直接除以分辨率了，可以换成"T y = (T)iY / ( parameters.getResolution() * parameters.getLy());"，这样和底下的IncomprFlowParam联系起来
    //Here y was directly divided by resolution, the line can be replaced by "T y = (T)iY / ( parameters.getResolution() * parameters.getLy());" 
    //to get the relation with IncomprFlowParam's parameter of ly, which means the height of the domain.
    return 4.*parameters.getLatticeU() * (y-y*y); //T u = parameters.getLatticeU()
}

/// A functional, used to initialize the velocity for the boundary conditions
//通过下面这个类，把poiseuilleVelocity和IncomprFlowParam结合起来,像这种类我们以后还会见到很多
//除了使用operator，我们还可以用virtual bool operator
/*By the class used below, the poiseuilleVelocity is linked with IncomprFlowParam, in the future we will see more.
*Additionally, we can also use virtual bool operator here for defining.
*/
template<typename T>
class PoiseuilleVelocity {
public:
    PoiseuilleVelocity(IncomprFlowParam<T> parameters_)
        : parameters(parameters_)
    { }
    /// This version of the operator returns the velocity only,
    ///    to instantiate the boundary condition.
    void operator()(plint iX, plint iY, Array<T,2>& u) const {
        u[0] = poiseuilleVelocity(iY, parameters);
        u[1] = T();
    }
    /// This version of the operator returns also a constant value for
    ///    the density, to create the initial condition.
    void operator()(plint iX, plint iY, T& rho, Array<T,2>& u) const {
        u[0] = poiseuilleVelocity(iY, parameters);
        u[1] = T();
        rho  = (T)1;
    }
private:
    IncomprFlowParam<T> parameters;
};
//以下这四行可以当作模板形式记忆，在随后的大括号里，描绘流域的所有信息
//Below 4 lines can be remembered as a template, in the subsequent brace written the whole domain's information.
void channelSetup (
        MultiBlockLattice2D<T,DESCRIPTOR>& lattice,
        IncomprFlowParam<T> const& parameters,
        OnLatticeBoundaryCondition2D<T,DESCRIPTOR>& boundaryCondition )
{
    // Create Velocity boundary conditions.
    boundaryCondition.setVelocityConditionOnBlockBoundaries(lattice);

/*如果需要的话，我们也可以在此定义
//If it is needed, we can also define:
*const plint nx = parameters.getNx();
*const plint ny = parameters.getNy();
*boundaryCondition.setVelocityConditionOnBlockBoundaries (
*           lattice, Box2D(nx-1, nx, 1, ny-2), boundary::outflow );
*来设置一个出口边界
//to set a outflow boundary.
*/



    // Specify the boundary velocity.
    //为边界设置速度
    setBoundaryVelocity (
            lattice, lattice.getBoundingBox(),
            PoiseuilleVelocity<T>(parameters) ); 

    // Create the initial condition.
    //此处尽管是为全部流场设置速度，但这儿只对之前定义为Dirichlet边界条件的格点生效
    //though here the velocity was attributed to all domain, only takes effect on nodes of Dirichlet boundary contion.
    initializeAtEquilibrium (
           lattice, lattice.getBoundingBox(), PoiseuilleVelocity<T>(parameters) );

    lattice.initialize();
}

void writeGifs(MultiBlockLattice2D<T,DESCRIPTOR>& lattice, plint iter)
{
    const plint imSize = 600;
    ImageWriter<T> imageWriter("leeloo");
    imageWriter.writeScaledGif(createFileName("u", iter, 6),
                               *computeVelocityNorm(lattice),
                               imSize, imSize );
}

int main(int argc, char* argv[]) {
    plbInit(&argc, &argv);

    global::directories().setOutputDir("./tmp/");

    // Use the class IncomprFlowParam to convert from
    //   dimensionless variables to lattice units, in the
    //   context of incompressible flows.
    IncomprFlowParam<T> parameters(
            (T) 1e-2,  // Reference velocity (the maximum velocity
                       //   in the Poiseuille profile) in lattice units.
            (T) 100.,  // Reynolds number//它的松弛参数会根据雷诺数自动计算
            //The relaxation time can be automatically calculated by Reynolds number here.
            100,       // Resolution of the reference length (channel height).
            2.,        // lx Channel length in dimensionless variables
            1.         // ly Channel height in dimensionless variables
    );
    const T imSave   = (T)0.1;  // Time intervals at which to save GIF
                                //   images, in dimensionless time units.
    const T maxT     = (T)3.1;  // Total simulation time, in dimensionless
                                //   time units.
    writeLogFile(parameters, "Poiseuille flow");

    MultiBlockLattice2D<T, DESCRIPTOR> lattice (
             parameters.getNx(), parameters.getNy(),
             new BGKdynamics<T,DESCRIPTOR>(parameters.getOmega()) );

//生成一个完全局部的矩形边界条件
//Create an entire local rectangular boundary condition.
    OnLatticeBoundaryCondition2D<T,DESCRIPTOR>*
        boundaryCondition = createLocalBoundaryCondition2D<T,DESCRIPTOR>();

    channelSetup(lattice, parameters, *boundaryCondition);

    // Main loop over time iterations.
    for (plint iT=0; iT*parameters.getDeltaT()<maxT; ++iT) {
        if (iT%parameters.nStep(imSave)==0 && iT>0) {
            pcout << "Saving Gif at time step " << iT << endl;
            writeGifs(lattice, iT);
        }
        // Execute lattice Boltzmann iteration.
        lattice.collideAndStream();
    }

    delete boundaryCondition;
}
