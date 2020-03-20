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

//#include "poiseuille.h"
//#include "poiseuille.hh"

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace plb;
using namespace std;

typedef double T;
#define DESCRIPTOR descriptors::D2Q9Descriptor



/// Velocity on the parabolic Poiseuille profile
T poiseuilleVelocity(plint iY, IncomprFlowParam<T> const& parameters) {
    T y = (T)iY / parameters.getResolution();
    return 4.*parameters.getLatticeU() * (y-y*y);
}
/// A functional, used to initialize the velocity for the boundary conditions
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

void channelSetup (
        MultiBlockLattice2D<T,DESCRIPTOR>& lattice,
        IncomprFlowParam<T> const& parameters,
        OnLatticeBoundaryCondition2D<T,DESCRIPTOR>& boundaryCondition )
{
    // Create Velocity boundary conditions.
    boundaryCondition.setVelocityConditionOnBlockBoundaries(lattice);

    // Specify the boundary velocity.
    setBoundaryVelocity (
            lattice, lattice.getBoundingBox(),
            PoiseuilleVelocity<T>(parameters) );

    // Create the initial condition.
    initializeAtEquilibrium (
           lattice, lattice.getBoundingBox(), PoiseuilleVelocity<T>(parameters) );

    lattice.initialize();
}





int main(int argc, char* argv[]) {
    plbInit(&argc, &argv);
    global::directories().setOutputDir("./tmp/");
    // Define numeric parameters.
    IncomprFlowParam<T> parameters (
        (T) 1e-2,  // uMax
        (T) 300.,  // Re
        299,        // N
        3.,        // lx
        1.         // ly 
    );

    writeLogFile(parameters, "Poiseuille flow");

    MultiBlockLattice2D<T, DESCRIPTOR> lattice (
              parameters.getNx(), parameters.getNy(),
              new BGKdynamics<T,DESCRIPTOR>(parameters.getOmega()) );

    OnLatticeBoundaryCondition2D<T,DESCRIPTOR>*
        boundaryCondition = createLocalBoundaryCondition2D<T,DESCRIPTOR>();

 MultiScalarField2D<bool> boolMask(parameters.getNx(), parameters.getNy());
    plb_ifstream ifile("1.dat");
    ifile >> boolMask;
    defineDynamics(lattice, boolMask, new BounceBack<T,DESCRIPTOR>, true);

   //createPoiseuilleBoundaries(lattice, parameters, *boundaryCondition);
   //lattice.initialize();
 channelSetup(lattice, parameters, *boundaryCondition);
    // Main loop over time iterations.
    for (plint iT=0; iT<40000; ++iT) {
        if (iT%1000==0) {
            pcout << "Writing image at dimensionless time " << iT*parameters.getDeltaT() << endl;
            ImageWriter<T> imageWriter("leeloo");
            imageWriter.writeScaledGif (
                    createFileName("velocity", iT, 6),
                    *computeVelocityNorm(lattice) );
            pcout << computeAverageEnergy(lattice) << endl;
        }

        // Lattice Boltzmann iteration step.
        lattice.collideAndStream();
    }

    delete boundaryCondition;
}
