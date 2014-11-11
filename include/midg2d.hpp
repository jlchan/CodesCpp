#ifndef __MIDG2D
#define __MIDG2D

#include <iostream>
#include <ostream>
#include <fstream>

/// High order elemental node coordinates and operators

#include "midg1d.hpp"
#include "mesh2d.hpp"

class midg2d;

class midg2d : public mesh2d {
public:

  /// Number of elements in total mesh
  int K;

  // Number of nodes on an edge, face, element
  int N,Nq,Nfp,Np;

  /// Number of cubature nodes in an element
  int Ncub, Ngauss, Ngauss_Nfaces;

  /// Nodal (W&B) derivative matrices
  fmatrix Dr, Ds;

  /// Nodal (cubature) derivative matrices
  fmatrix cDr, cDs;

  /// Weak cubature derivative matrices and interpolation matrix
  fmatrix cMinvDrTW, cMinvDsTW, cMinvIcTW, cinterp;

  /// Gauss interpolation lift matrix
  fmatrix ginterp, gMinvIcTW;

  /// vertex interpolation matrix
  fmatrix vinterp;

  /// Physical coordinates of W&B nodes
  fmatrix x,y;

  /// Physical coordinates of cubature nodes
  fmatrix cx,cy;

  /// Physical coordinates of element surface Gauss nodes
  fmatrix gx,gy;

  /// Element local coordinates of Warp & Blend Nodes
  fmatrix r,s;

  /// Element local coordinates of cubature nodes
  fmatrix cr,cs,cw;

  /// Element local coordinates of Gauss nodes
  fmatrix gr,gs,gw;

  ///
  fmatrix V, Vinv;

  /// Geometric factors
  fmatrix rx, sx, ry,sy, J;

  /// Surface geometry info
  fmatrix nx, ny, sJ, Fscale;

  // Constructor: calls base class to load mesh in
  midg2d(setupAide &setup);

  midg2d();

  void startup();
  /**
   * This function sets up the W&B nodes (local coordinates) and the
   * collocation differentiation matrix.
   *
   */
  void initNodes(int inN, int inC, int inG, int inP);

  /**
   * This function finds the W&B nodes in physical space and then computes geometric factors.
   *
   */
  void geometry();

};

//  datafloat factorial(int n);

datafloat Warpfactor(const int p, fmatrix &rout);

void Nodes2D(const int p, fmatrix &x, fmatrix &y);

void xytors(const fmatrix &x, const fmatrix &y,
	    fmatrix &r, fmatrix &s);

datafloat Simplex2DP(datafloat aout, datafloat bout,
		     int i, int j);

void GradSimplex2DP(datafloat aout, datafloat bout,
		    int id, int jd,
		    datafloat &dmodedr,
		    datafloat &dmodeds);

void rstoab(const fmatrix &r, const fmatrix &s,
	    fmatrix &a, fmatrix &b);

void Vandermonde2D(int p, fmatrix &rout,
		   fmatrix &sout,
		   fmatrix &Vout);

void GradVandermonde2D(int p, fmatrix &rout,
		       fmatrix &sout,
		       fmatrix &Vrout,
		       fmatrix &Vsout);

void basis2D(int p, fmatrix &rout, fmatrix &sout,
	     fmatrix &Vout, fmatrix &Vrout,
	     fmatrix &Vsout);

void Dmatrices2D(const int p, const fmatrix &r,
		 const fmatrix &s,
		 const fmatrix &V,
		 fmatrix &Dr,
		 fmatrix &Ds);

void Lift2D(fmatrix &Lift);


#endif
