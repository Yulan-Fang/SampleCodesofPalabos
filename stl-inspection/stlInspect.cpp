/* This file is part of the Palabos library.
 *
 * Copyright (C) 2011-2017 FlowKit Sarl
 * Route d'Oron 2
 * 1010 Lausanne, Switzerland
 * E-mail contact: contact@flowkit.com
 *
 * The most recent release of Palabos can be downloaded at
 * <http://www.palabos.org/>
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

/** \file
  * Flow in a lid-driven 3D cavity. The cavity is square and has no-slip walls,
  * except for the top wall which is diagonally driven with a constant
  * velocity. The benchmark is challenging because of the velocity
  * discontinuities on corner nodes.
  **/

#include "palabos3D.h"
#include "palabos3D.hh"
#include "palabos2D.h"
#include "palabos2D.hh"
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>

using namespace plb;
using namespace std;

typedef double T;
#define DESCRIPTOR descriptors::D3Q19Descriptor

void cavitySetup( MultiBlockLattice3D<T,DESCRIPTOR>& lattice,
                  IncomprFlowParam<T> const& parameters,
                  OnLatticeBoundaryCondition3D<T,DESCRIPTOR>& boundaryCondition )
{
    const plint nx = parameters.getNx();
    const plint ny = parameters.getNy();
    const plint nz = parameters.getNz();
    Box3D topLid = Box3D(0, nx-1, ny-1, ny-1, 0, nz-1);
    Box3D everythingButTopLid = Box3D(0, nx-1, 0, ny-2, 0, nz-1);

    boundaryCondition.setVelocityConditionOnBlockBoundaries(lattice, lattice.getBoundingBox(), boundary::dirichlet);

    T u = std::sqrt((T)2)/(T)2 * parameters.getLatticeU();
    initializeAtEquilibrium(lattice, everythingButTopLid, (T) 1., Array<T,3>((T)0.,(T)0.,(T)0.) );
    initializeAtEquilibrium(lattice, topLid, (T) 1., Array<T,3>(u,(T)0.,u) );
    //setBoundaryVelocity(lattice, topLid, Array<T,3>(u,(T)0.,u) );
    setBoundaryVelocity(lattice, topLid, Array<T,3>(u,(T)0.,u/2.0) );

    lattice.initialize();
}

template<class BlockLatticeT>
void writeGifs(BlockLatticeT& lattice,
               IncomprFlowParam<T> const& parameters, plint iter)
{
    const plint imSize = 600;
    const plint nx = parameters.getNx();
    const plint ny = parameters.getNy();
    const plint nz = parameters.getNz();
    const plint zComponent = 2;

    Box3D slice(0, nx-1, 0, ny-1, nz/2, nz/2);
    ImageWriter<T> imageWriter("leeloo");

    imageWriter.writeScaledGif( createFileName("uz", iter, 6),
                                *computeVelocityComponent (lattice, slice, zComponent),
                                imSize, imSize );

    imageWriter.writeScaledGif( createFileName("uNorm", iter, 6),
                                *computeVelocityNorm (lattice, slice),
                                imSize, imSize );
    imageWriter.writeScaledGif( createFileName("omega", iter, 6),
                                *computeNorm(*computeVorticity (
                                        *computeVelocity(lattice) ), slice ),
                                imSize, imSize );
}

template<class BlockLatticeT>
void writeVTK(BlockLatticeT& lattice,
              IncomprFlowParam<T> const& parameters, plint iter)
{
    T dx = parameters.getDeltaX();
    T dt = parameters.getDeltaT();

    //VtkImageOutput3D<T> vtkOut(createFileName("vtk", iter, 6), dx);
    ParallelVtkImageOutput3D<T> vtkOut(createFileName("vtk", iter, 6), 3, dx);

    vtkOut.writeData<3,float>(*computeVelocity(lattice), "velocity", dx/dt);
    vtkOut.writeData<float>(*computeVelocityNorm(lattice), "velocityNorm", dx/dt);
    vtkOut.writeData<3,float>(*computeVorticity(*computeVelocity(lattice)), "vorticity", 1./dt);
}

void writeSurface(TriangleSet<T> triangleSet, plint iT)
{
    triangleSet.writeBinarySTL(createFileName("surface_", iT, 8) + ".stl");
}

int main(int argc, char* argv[]) {

    plbInit(&argc, &argv);
    global::directories().setOutputDir("./tmp/");
    //defaultMultiBlockPolicy3D().toggleBlockingCommunication(true);

    IncomprFlowParam<T> parameters(
            (T) 1e-2,  // uMax
            (T) 10.,   // Re
            1,        // N
            100.,        // lx
            100.,        // ly
            100.         // lz
    );
    const T logT     = (T)1/(T)10;
    const T imSave   = (T)1/(T)10;
    const T vtkSave  = (T)1/(T)100;
    const T maxT     = (T)1.1;

    pcout << "omega= " << parameters.getOmega() << std::endl;
    writeLogFile(parameters, "3D diagonal cavity");

    T omega = parameters.getOmega();

    MultiBlockLattice3D<T, DESCRIPTOR> lattice (
            parameters.getNx(), parameters.getNy(), parameters.getNz(),
            new BGKdynamics<T,DESCRIPTOR>(omega) );

//% Inspection of the surface_
    std::vector<Array<T,3> > vertices;
    std::vector<T> areas;
    pcout << "reading the stl file." << std::endl;
    std::string meshFileName = "input.stl"; //The target stl file
    TriangleSet<T> triangleSet(meshFileName, DBL);

    Cuboid<T> bCuboid = triangleSet.getBoundingCuboid();
    Array<T,3> objectCenter = (T) 0.5 * (bCuboid.lowerLeftCorner + bCuboid.upperRightCorner);
    triangleSet.translate(-objectCenter);
    pcout << "The mesh center was moved to the (0,0,0)." << std::endl;
    triangleSet.scale((T) 1 / 100); // convert to LB unit
    triangleSet.scale((T) 100 / 100);//scaled to 1x%
    pcout << "The mesh scaled" << std::endl;
    DEFscaledMesh<T> *fishDef = new DEFscaledMesh<T>(triangleSet, 0, 0, 0, Dot3D(0, 0, 0));
    pcout << "The rectangle has " << fishDef->getMesh().getNumVertices() << " vertices and " <<
    fishDef->getMesh().getNumTriangles() << " triangles." << std::endl;
    for (pluint iVertex = 0; iVertex < (pluint) fishDef->getMesh().getNumVertices(); iVertex++) {
      vertices.push_back(fishDef->getMesh().getVertex(iVertex));
      areas.push_back(fishDef->getMesh().computeVertexArea(iVertex));
}
delete fishDef;
// It should be noticed that no FSI algorithm implemented here,
// so the stl file won't affect the flow


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    OnLatticeBoundaryCondition3D<T,DESCRIPTOR>* boundaryCondition
        //= createInterpBoundaryCondition3D<T,DESCRIPTOR>();
        = createLocalBoundaryCondition3D<T,DESCRIPTOR>();

    cavitySetup(lattice, parameters, *boundaryCondition);

    T previousIterationTime = T();
    // Loop over main time iteration.
    writeSurface(triangleSet, 0);
    writeVTK(lattice, parameters, 0);
    for (plint iT=0; iT<parameters.nStep(maxT); ++iT) {
        global::timer("mainLoop").restart();
        if (iT%parameters.nStep(vtkSave)==0 && iT>0) {
            pcout << "Saving VTK file and stl file..." << endl;
            writeVTK(lattice, parameters, iT);
            writeSurface(triangleSet, iT);
        }

        if (iT%parameters.nStep(logT)==0) {
            pcout << "step " << iT
                  << "; t=" << iT*parameters.getDeltaT();
        }

        // Execute a time iteration.
        lattice.collideAndStream();

        // Access averages from internal statistics ( their value is defined
        //   only after the call to lattice.collideAndStream() )
        if (iT%parameters.nStep(logT)==0) {
            pcout << "; av energy="
                  << setprecision(10) << getStoredAverageEnergy<T>(lattice)
                  << "; av rho="
                  << setprecision(10) << getStoredAverageDensity<T>(lattice) << endl;
            pcout << "Time spent during previous iteration: "
                  << previousIterationTime << endl;
        }

        previousIterationTime = global::timer("mainLoop").stop();
    }

    delete boundaryCondition;
}
