#ifndef __MESH
#define __MESH

/// source of the exit function
// #include "headers2d.hpp"

#include "matrix.hpp"
#include "setupAide.hpp"

class mesh2d;

enum BCnums {
  BC_Wall = 1,
  BC_Inflow = 2,
  BC_Outflow = 3,
  BC_Radiating = 4,
  BC_Dirichlet = 5
};


class mesh2d{

protected:

  int Nverts;
  int Nfaces;
  int Nvertsperface;

  /// original ordering of nodes in gmsh file
  matrix <datafloat> VX;   /// node x-coordinates: row vector for V
  matrix <datafloat> VY;   /// node y-coordinates: row vector for VY

  matrix <datafloat> VB; /// node bathymetry

  matrix <int>    EToV; /// element to vertex connectivity array

  matrix <int>    EToE; /// element to element connectivity array
  matrix <int>    EToF; /// element to element connectivity array

  matrix <int>    BCType; /// list of boundary types

  int nboun;
  matrix <int> BSID;
public:

  mesh2d();

  /// constructor: populates mesh with stuff from file named filename
  mesh2d(setupAide &setup);

  void connect();

  void diagnose();

  void debug(const char *name);

};

#endif
