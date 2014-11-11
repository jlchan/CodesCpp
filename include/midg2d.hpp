#ifndef __MIDG2D
#define __MIDG2D

#include <iostream>
#include <ostream>
#include <fstream>

//#include "headers2d.hpp"
//#define DEBUG_MIDG
/// High order elemental node coordinates and operators

#include "mesh2d.hpp"

class midg2d;

class midg2d : public mesh2d {
public:

  /// Number of elements in total mesh
  int K;

  int padding_factor;

  int thread_padding_factor;

  // Number of nodes on an edge, face, element
  int N,Nq,Nfp,Np,Npad;

  /// Number of cubature nodes in an element
  int Ncub, NcubPad;

  /// Number of Gauss integration nodes on a face
  int Ngauss, Ngauss_Nfaces, Ngauss_NfacesPad;

  /// Number of triangles for plotting
  int Nptris, Npp;
  imatrix ptris;
  fmatrix pr, ps, pinterp, px, py;

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

  /// interpolation matrices for positivity preservation
  fmatrix checkinterp, V1, gV1, P1, P0;

  ///
  fmatrix V, Vinv;

  fmatrix V1Nodal, P1Nodal, gV1Nodal;

  /// Geometric factors
  fmatrix rx, sx, ry,sy, J;

  /// Surface geometry info
  fmatrix nx, ny, sJ, Fscale;

  /// Number of levels in mrab
  int Nlevels;

  /// level information
  imatrix levflag;

  /// list of list of elements at each level
  matrix < imatrix > ks;

  /// list of coarse neighbors at each level
  matrix < imatrix > kcoarse;

  /// list of ?
  matrix < imatrix > kids;

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

  //  datafloat factorial(int n);

  datafloat JacobiP(datafloat xout, datafloat alpha, datafloat beta, int p);

  datafloat GradJacobiP(datafloat xout, datafloat alpha, datafloat beta,  int p);

  datafloat SimplexP(datafloat aout, datafloat bout, int i, int j);

  void GradSimplexP(datafloat aout, datafloat bout,
		    int id, int jd,
		    datafloat &dmodedr, datafloat &dmodeds);

  void rstoab(fmatrix &r, fmatrix &s, fmatrix &a, fmatrix &b);

  void Vandermonde(int p, fmatrix &rout, fmatrix &sout, fmatrix &Vout);

  void GradVandermonde(int p, fmatrix &rout, fmatrix &sout,
		       fmatrix &Vrout, fmatrix &Vsout);

  void basis(int p, fmatrix &rout, fmatrix &sout, fmatrix &Vout,
	     fmatrix &Vrout, fmatrix &Vsout);

  datafloat build_levels(fmatrix &all_dt, int maxNlevels);
};

#endif
