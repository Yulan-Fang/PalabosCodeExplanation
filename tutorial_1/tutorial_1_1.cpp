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

/* Code 1.1 in the Palabos tutorial
 */

#include "palabos2D.h"
#include "palabos2D.hh"
#include <iostream>
#include <iomanip>

using namespace plb;
using namespace std;

//T是双精度浮点数类型
//double-precision floating type T
typedef double T;

//使用D2Q9格子表示器
//Use D2Q9 Lattice descriptor
#define DESCRIPTOR plb::descriptors::D2Q9Descriptor

// Initialize the lattice at zero velocity and constant density, except
//   for a slight density excess on a square sub-domain.
// 以常数密度和零速度定义lattice，而一小块子域的密度要稍稍高一些
void defineInitialDensityAtCenter(MultiBlockLattice2D<T,DESCRIPTOR>& lattice)
{
    // The lattice is of size nx-by-ny
// 用plint定义常量nx和ny，分别为流域的x和y长度
// Use "plint" to define constant nx and ny, represent the length of x and y of the domain, respectively.
    const plint nx = lattice.getNx();
    const plint ny = lattice.getNy();

    // Create a Box2D which describes the location of cells with a slightly
    //  higher density.
// 下面使用plint定义三个变量
//use "plint" to define three variables below
//1）方形边长的一半
//1）The half of the length of the square 
    plint centralSquareRadius = nx/6; 
    
//2）流域x轴三分之一的位置，此处意为子域方形的中心点x轴坐标
//2）At the position of 1/3 of x axis, means the x coordinate of square center
    plint centerX = nx/3;
//3) plint the y coordinate of the center
    plint centerY = ny/4;

//使用Box2D来定义centralSquare，需输入x1,x2,y1,y2
//其中x轴从x1至x2，y轴从y1到y2，四个点确定2D域
//此处x1 = centerX - centralSquareRadius，x2 = centerX + centralSquareRadius
//Use Box2D to define centralSquare, as it needs input of x1,x2,y1,y2. 
//In x axis from x1 to x2, and in y axis from y1 to y2, then 4 points determine a 2D domain
//here x1 = centerX - centralSquareRadius，x2 = centerX + centralSquareRadius
    Box2D centralSquare (
            centerX - centralSquareRadius, centerX + centralSquareRadius,
            centerY - centralSquareRadius, centerY + centralSquareRadius );

    // All cells have initially density rho ...
// 所有的cell初始化密度rho
    T rho0 = 1.;
    // .. except for those in the box "centralSquare" which have density
    //    rho+deltaRho
//为流域中部的方形子域准备一个稍微大一点的密度
    T deltaRho = 1.e-4;

//像规模较大的数据，如LB模拟的粒子群，都储存在palabos的BlockLattice类型里
//像Array<T,nDim>一般是小型数据，此处为T数据类型（在程序最前已定义为双精度浮点型）
//nDim为尺寸，此处2表示这个固定尺寸的array有两个元素
//你也可以在定义这2或3元素Array的同时，储存一些物理向量进去
//For large size data, like particles in Lattice Boltzmann simulation, were stored in BlockLattice type of Palabos
//As Array<T,nDim> generally small, here it is of T data type(In the beginning of the code T was difined as double-precision floating point scalars)
//nDim is the size, here is 2, which represents this fixed-size array has two elements
// You can also store small physical vectors, such as velocity, when you define the array
    Array<T,2> u0((T)0,(T)0);

    // Initialize constant density everywhere.
//初始化常数密度，initializeAtEquilibrium(lattice, domain, rho, velocity)，在该区域的cells中初始化平衡分布
//initializeAtEquilibrium(lattice, domain, rho, velocity), initialize constant density rho to the cells in the domain
    initializeAtEquilibrium (
           lattice, lattice.getBoundingBox(), rho0, u0 );

    // And slightly higher density in the central box.
//中间方形子域的密度稍微高一些
    initializeAtEquilibrium (
           lattice, centralSquare, rho0 + deltaRho, u0 );
//lattice初始化
    lattice.initialize();
}//Here "void defineInitialDensityAtCenter(MultiBlockLattice2D<T,DESCRIPTOR>& lattice)" part is finished
//记住这个defineInitialDensityAtCenter，后面还会出现
//Remember this "defineInitialDensityAtCenter", later will reappear


//下面这行代码，一般通过argc和argv在cpp中用于访问程序参数
//在程序的任何位置都可以访问这两个全局对象
//Usually, we access to parameters in CPP programs through argc and argv
//In anywhere of the program can we access these two global parameters
int main(int argc, char* argv[]) {
    plbInit(&argc, &argv);

//设定文件输出位置，即打开tmp文件夹
//output file directory, that is to open floder tmp
    global::directories().setOutputDir("./tmp/");

    const plint maxIter = 1000; // Iterate during 1000 steps.迭代次数
    const plint nx = 600;       // Choice of lattice dimensions.流域x轴方向lattice格子数
    const plint ny = 600;
    const T omega = 1.;        // Choice of the relaxation parameter松弛时间

//初始化一个nx,ny，BGK模型的T数据类型（即双精度浮点型）的D2Q9（DESCRIPTOR在程序最上面已定义）的lattice
//initialize a nx,ny size, BGK dynamics, type T(Double-precision floating), D2Q9 lattice
//BGKdynamics: The well known BGK single-relaxation-time model, which works on D2Q9, D3Q15, D3Q19 and D3Q27 ——from Palabos User Guide
    MultiBlockLattice2D<T, DESCRIPTOR> lattice (
           nx, ny, new BGKdynamics<T,DESCRIPTOR>(omega) );

//true或false来调整循环边界
//Select true or false to turn on or off the periodic boundaries
    lattice.periodicity().toggleAll(true); // Use periodic boundaries.

//还记得这个defineInitialDensityAtCenter吗？
//You remember this defineInitialDensityAtCenter?
    defineInitialDensityAtCenter(lattice);

    // Main loop over time iterations.
// 主要的迭代过程
    for (plint iT=0; iT<maxIter; ++iT) {
        if (iT%40==0) {  // Write an image every 40th time step. 以40取余为0时
            pcout << "Writing GIF file at iT=" << iT << endl;
            // Instantiate an image writer with the color map "leeloo".
            ImageWriter<T> imageWriter("leeloo");
// 初始化绘图，颜色表leeloo
// leeloo可以换成earth, water, air, fire
            // Write a GIF file with colors rescaled to the range of values
            //   in the matrix
// 启动绘图
            imageWriter.writeScaledGif (
                    createFileName("u", iT, 6),
                    *computeVelocityNorm(lattice) );
        }//绘图结束
        // Execute lattice Boltzmann iteration.
// 执行碰撞和流动
        lattice.collideAndStream();
    }//下一轮循环，直到iT达到maxIter
}

// Code explanation by Yulan Fang 
// Error correction please send to ahdhfang@hotmail.com
// ——March,9 2020 at Siwa, Egypt.