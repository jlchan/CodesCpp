#include <stdlib.h>
#include <ostream>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>
#include <assert.h>
#include "midg2d.hpp"

/// High order elemental node coordinates and operators

// Constructor: calls base class to load mesh in
midg2d::midg2d(setupAide &setup) : mesh2d(setup){

  /// number of nodes on each edge, face, element
  if(!setup.getArgs("POLYNOMIAL DEGREE", N)) throw -1;

  Nq  = N+1;
  Nfp = (N+1);
  Np  = (N+1)*(N+2)/2;

  K = EToV.nrows();

  connect();

  initNodes(N, 2*N, 2*N, 2*N);

  geometry();
}

midg2d::midg2d(){
  N = 0;
  Nq = 0;
  Nfp = 0;
  Np = 0;
}

void midg2d::startup(){

  /// build element-element connectivity information
  connect();

  /// set up physical W&B node locations
  initNodes(N, 2*N, 2*N, 2*N);

  /// compute geometric factors and surface normals
  geometry();
}

/**
 * This function sets up the W&B nodes (local coordinates) and the
 * collocation differentiation matrix.
 *
 */
void midg2d::initNodes(int inN, int inC, int inG, int inP){

  int N = inN, G = inG, P = inP, C = inC;

#include "data.h"

  /// Build coordinates
  Ncub   = ref_Ncub[inC];
  Ngauss = ref_Ngauss[inG];

  std::cout << " Np = " << Np << std::endl;
  std::cout << " Ngauss_Nfaces = " << Ngauss_Nfaces << std::endl;
  std::cout << " Ncub = " << Ncub << std::endl;

#ifdef DEBUG_MIDG
  std::cout << " Ncub = " << Ncub << std::endl;
  std::cout << " Ngauss = " << Ngauss << std::endl;
#endif

  /// resize Gauss node matrices and data
  gr.resize(Ngauss*Nfaces, 1);
  gs.resize(Ngauss*Nfaces, 1);
  gw.resize(Ngauss*Nfaces, 1);

  /// resize cubature matrices and data
  cr.resize(Ncub, 1);
  cs.resize(Ncub, 1);
  cw.resize(Ncub, 1);

  /// resize data for Warp&Blend nodes
  r.resize(Np, 1);
  s.resize(Np, 1);

  double *ref_r, *ref_s;
  switch(N){
  case 1: ref_r = ref_r_01; ref_s = ref_s_01; break;
  case 2: ref_r = ref_r_02; ref_s = ref_s_02; break;
  case 3: ref_r = ref_r_03; ref_s = ref_s_03; break;
  case 4: ref_r = ref_r_04; ref_s = ref_s_04; break;
  case 5: ref_r = ref_r_05; ref_s = ref_s_05; break;
  case 6: ref_r = ref_r_06; ref_s = ref_s_06; break;
  case 7: ref_r = ref_r_07; ref_s = ref_s_07; break;
  case 8: ref_r = ref_r_08; ref_s = ref_s_08; break;
  }

  /// element local coordinates of Warp&Blend nodes
  for(int i=1;i<=Np;++i){
    r(i) = ref_r[i-1];
    s(i) = ref_s[i-1];
  }

  double *ref_cr, *ref_cs, *ref_cw;
  switch(C){
  case  1: ref_cr = ref_cr_01; ref_cs = ref_cs_01; ref_cw = ref_cw_01; break;
  case  2: ref_cr = ref_cr_02; ref_cs = ref_cs_02; ref_cw = ref_cw_02; break;
  case  3: ref_cr = ref_cr_03; ref_cs = ref_cs_03; ref_cw = ref_cw_03; break;
  case  4: ref_cr = ref_cr_04; ref_cs = ref_cs_04; ref_cw = ref_cw_04; break;
  case  5: ref_cr = ref_cr_05; ref_cs = ref_cs_05; ref_cw = ref_cw_05; break;
  case  6: ref_cr = ref_cr_06; ref_cs = ref_cs_06; ref_cw = ref_cw_06; break;
  case  7: ref_cr = ref_cr_07; ref_cs = ref_cs_07; ref_cw = ref_cw_07; break;
  case  8: ref_cr = ref_cr_08; ref_cs = ref_cs_08; ref_cw = ref_cw_08; break;
  case  9: ref_cr = ref_cr_09; ref_cs = ref_cs_09; ref_cw = ref_cw_09; break;
  case 10: ref_cr = ref_cr_10; ref_cs = ref_cs_10; ref_cw = ref_cw_10; break;
  case 11: ref_cr = ref_cr_11; ref_cs = ref_cs_11; ref_cw = ref_cw_11; break;
  case 12: ref_cr = ref_cr_12; ref_cs = ref_cs_12; ref_cw = ref_cw_12; break;
  case 13: ref_cr = ref_cr_13; ref_cs = ref_cs_13; ref_cw = ref_cw_13; break;
  case 14: ref_cr = ref_cr_14; ref_cs = ref_cs_14; ref_cw = ref_cw_14; break;
  case 15: ref_cr = ref_cr_15; ref_cs = ref_cs_15; ref_cw = ref_cw_15; break;
  case 16: ref_cr = ref_cr_16; ref_cs = ref_cs_16; ref_cw = ref_cw_16; break;
  case 17: ref_cr = ref_cr_17; ref_cs = ref_cs_17; ref_cw = ref_cw_17; break;
  case 18: ref_cr = ref_cr_18; ref_cs = ref_cs_18; ref_cw = ref_cw_18; break;
  case 19: ref_cr = ref_cr_19; ref_cs = ref_cs_19; ref_cw = ref_cw_19; break;
  }

  /// element local coordinates of cubature nodes
  for(int i=1;i<=Ncub;++i){
    cr(i) = ref_cr[i-1];
    cs(i) = ref_cs[i-1];
    cw(i) = ref_cw[i-1];
  }


  double *ref_gr, *ref_gs, *ref_gw;
  switch(G){
  case  1: ref_gr = ref_gr_01; ref_gs = ref_gs_01; ref_gw = ref_gw_01; break;
  case  2: ref_gr = ref_gr_02; ref_gs = ref_gs_02; ref_gw = ref_gw_02; break;
  case  3: ref_gr = ref_gr_03; ref_gs = ref_gs_03; ref_gw = ref_gw_03; break;
  case  4: ref_gr = ref_gr_04; ref_gs = ref_gs_04; ref_gw = ref_gw_04; break;
  case  5: ref_gr = ref_gr_05; ref_gs = ref_gs_05; ref_gw = ref_gw_05; break;
  case  6: ref_gr = ref_gr_06; ref_gs = ref_gs_06; ref_gw = ref_gw_06; break;
  case  7: ref_gr = ref_gr_07; ref_gs = ref_gs_07; ref_gw = ref_gw_07; break;
  case  8: ref_gr = ref_gr_08; ref_gs = ref_gs_08; ref_gw = ref_gw_08; break;
  case  9: ref_gr = ref_gr_09; ref_gs = ref_gs_09; ref_gw = ref_gw_09; break;
  case 10: ref_gr = ref_gr_10; ref_gs = ref_gs_10; ref_gw = ref_gw_10; break;
  case 11: ref_gr = ref_gr_11; ref_gs = ref_gs_11; ref_gw = ref_gw_11; break;
  case 12: ref_gr = ref_gr_12; ref_gs = ref_gs_12; ref_gw = ref_gw_12; break;
  case 13: ref_gr = ref_gr_13; ref_gs = ref_gs_13; ref_gw = ref_gw_13; break;
  case 14: ref_gr = ref_gr_14; ref_gs = ref_gs_14; ref_gw = ref_gw_14; break;
  case 15: ref_gr = ref_gr_15; ref_gs = ref_gs_15; ref_gw = ref_gw_15; break;
  case 16: ref_gr = ref_gr_16; ref_gs = ref_gs_16; ref_gw = ref_gw_16; break;
  case 17: ref_gr = ref_gr_17; ref_gs = ref_gs_17; ref_gw = ref_gw_17; break;
  case 18: ref_gr = ref_gr_18; ref_gs = ref_gs_18; ref_gw = ref_gw_18; break;
  case 19: ref_gr = ref_gr_19; ref_gs = ref_gs_19; ref_gw = ref_gw_19; break;
  }

  /// element local coordinates of cubature nodes
  for(int i=1;i<=Ngauss*Nfaces;++i){
    gr(i) = ref_gr[i-1];
    gs(i) = ref_gs[i-1];
    gw(i) = ref_gw[i-1];
  }

#if 0
  std::cout << "r = " << r << std::endl;
  std::cout << "s = " << s << std::endl;

  std::cout << "cr = " << cr << std::endl;
  std::cout << "cs = " << cs << std::endl;
  std::cout << "cw = " << cw << std::endl;

  std::cout << "gr = " << gr << std::endl;
  std::cout << "gs = " << gs << std::endl;
  std::cout << "gw = " << gw << std::endl;
#endif
  /// test recurrence
  fmatrix  Vr,  Vs;
  fmatrix cV, cVr, cVs;
  fmatrix gV, gVr, gVs;
  fmatrix vr, vs, vV, vVr, vVs;

  vr.resize(3,1);
  vs.resize(3,1);

  for(int i=1; i<=3; ++i){
    vr(i) = ref_r_01[i-1];
    vs(i) = ref_r_01[i-1];
  }

  basis2D(1, vr, vs, vV, vVr, vVs);
  basis2D(N,  r,  s,  V,  Vr,  Vs);
  basis2D(N, cr, cs, cV, cVr, cVs);
  basis2D(N, gr, gs, gV, gVr, gVs);

  fmatrix MMinv = V*(V.transpose());

  fmatrix cW(Ncub,Ncub);
  for(int n=1;n<=Ncub;++n){
    cW(n,n) = cw(n);
  }

  fmatrix gW(Ngauss*Nfaces,Ngauss*Nfaces);
  for(int n=1;n<=Ngauss*Nfaces;++n){
    gW(n,n) = gw(n);
  }

  Vinv = V.inverse();

  cinterp = cV*Vinv;
  ginterp = gV*Vinv;
  vinterp = vV*Vinv;

  cMinvIcTW = MMinv*(cinterp.transpose()*cW);
  gMinvIcTW = MMinv*(ginterp.transpose()*gW);

  cDr = cVr*Vinv;
  cDs = cVs*Vinv;

  Dr = Vr*Vinv;
  Ds = Vs*Vinv;

  cMinvDrTW = MMinv*(cDr.transpose()*cW);
  cMinvDsTW = MMinv*(cDs.transpose()*cW);
}

/**
 * This function finds the W&B nodes in physical space and then computes geometric factors.
 *
 **/
void midg2d::geometry(){

  /// physical node coordinates
  x.resize(Np, K);
  y.resize(Np, K);

  /// cubature node coordinates
  cx.resize(Ncub, K);
  cy.resize(Ncub, K);

  /// element boundary Gauss node coordinates
  gx.resize(Ngauss*Nfaces, K);
  gy.resize(Ngauss*Nfaces, K);

  /// Volume geometric factors
  rx.resize(1, K); ry.resize(1, K);
  sx.resize(1, K); sy.resize(1, K);
  J.resize(1, K);

  /*! Surface geometric factors */
  nx.resize(Nfaces, K);
  ny.resize(Nfaces, K);
  sJ.resize(Nfaces, K);
  Fscale.resize(Nfaces, K);

  for(int k=1;k<=K;++k){

    double x1 = VX[EToV(k,1)];
    double x2 = VX[EToV(k,2)];
    double x3 = VX[EToV(k,3)];
    double y1 = VY[EToV(k,1)];
    double y2 = VY[EToV(k,2)];
    double y3 = VY[EToV(k,3)];

    // set W&B physical coordinates for element k
    for(int i=1;i<=Np;++i){
      x(i,k) = 0.5*( -(r(i)+s(i))*x1 + (1+r(i))*x2 + (1+s(i))*x3 );
      y(i,k) = 0.5*( -(r(i)+s(i))*y1 + (1+r(i))*y2 + (1+s(i))*y3 );
    }

    // set cubature physical coordinates for element k
    for(int i=1;i<=Ncub;++i){
      cx(i,k) = 0.5*( -(cr(i)+cs(i))*x1 + (1+cr(i))*x2 + (1+cs(i))*x3 );
      cy(i,k) = 0.5*( -(cr(i)+cs(i))*y1 + (1+cr(i))*y2 + (1+cs(i))*y3 );
    }

    // set Gauss face nodes physical coordinates for element k
    for(int i=1;i<=Ngauss*Nfaces;++i){
      gx(i,k) = 0.5*( -(gr(i)+gs(i))*x1 + (1+gr(i))*x2 + (1+gs(i))*x3 );
      gy(i,k) = 0.5*( -(gr(i)+gs(i))*y1 + (1+gr(i))*y2 + (1+gs(i))*y3 );
    }

    double xr = 0.5*(x2-x1);
    double xs = 0.5*(x3-x1);
    double yr = 0.5*(y2-y1);
    double ys = 0.5*(y3-y1);

    J(k) = xr*ys-xs*yr;

    rx(k) =  ys/J(k);
    sx(k) = -yr/J(k);
    ry(k) = -xs/J(k);
    sy(k) =  xr/J(k);

    nx(1,k) = yr;    ny(1,k) = -xr;
    nx(2,k) = ys-yr; ny(2,k) = -xs+xr;
    nx(3,k) =-ys;    ny(3,k) =  xs;
    for(int f=1;f<=Nfaces;++f){
      sJ(f,k) = sqrt(nx(f,k)*nx(f,k)+ny(f,k)*ny(f,k));
      nx(f,k) = nx(f,k)/sJ(f,k);
      ny(f,k) = ny(f,k)/sJ(f,k);
      Fscale(f,k) = sJ(f,k)/J(k);
    }
  }

  std::cout << " x range: [" << x.minentry() << "," << x.maxentry() << "]" << std::endl;
  std::cout << " y range: [" << y.minentry() << "," << y.maxentry() << "]" << std::endl;
  std::cout << " sJ range: [" << sJ.minentry() << "," << sJ.maxentry() << "]" << std::endl;
  std::cout << " Fscale range: [" << Fscale.minentry() << "," << Fscale.maxentry() << "]" << std::endl;

}

datafloat Warpfactor(const int p, fmatrix &rout){
}

void Nodes2D(const int p, fmatrix &xout,
	     fmatrix &yout){
}

void xytors(const fmatrix &x, const fmatrix &y,
	    fmatrix &r, fmatrix &s);

datafloat Simplex2DP(datafloat aout, datafloat bout,
		     int i, int j){

  /// SimplexP(a,b,i,j);
  /// Purpose : Evaluate 2D orthonormal polynomial
  ///           on simplex at (a,b) of order (i,j).

  datafloat h1 = JacobiP(aout,0,0,i);
  datafloat h2 = JacobiP(bout,2*i+1,0,j);
  datafloat P = sqrt(2.0)*h1*h2*pow(1-bout,i);
  return P;
}


void GradSimplex2DP(datafloat aout, datafloat bout,
		    int id, int jd,
		    datafloat &dmodedr,
		    datafloat &dmodeds){

  /// function [dmodedr, dmodeds] = GradSimplex2DP(a,b,id,jd)
  /// Purpose: Return the derivatives of the modal basis (id,jd)
  ///          on the  simplex at (a,b).

  datafloat  fa = JacobiP(aout, 0, 0, id);
  datafloat dfa = GradJacobiP(aout, 0, 0, id);
  datafloat  gb = JacobiP(bout, 2*id+1,0, jd);
  datafloat dgb = GradJacobiP(bout, 2*id+1,0, jd);

  /// r-derivative
  /// d/dr = da/dr d/da + db/dr d/db = (2/(1-s)) d/da = (2/(1-b)) d/da
  dmodedr = dfa*gb;
  if(id>0)
    dmodedr = dmodedr*pow(0.5*(1-bout),(id-1));

  /// s-derivative
  /// d/ds = ((1+a)/2)/((1-b)/2) d/da + d/db
  dmodeds = dfa*(gb*(0.5*(1+aout)));
  if(id>0)
    dmodeds = dmodeds*pow(0.5*(1-bout),(id-1));

  datafloat tmp = dgb*pow(0.5*(1-bout),id);
  if(id>0)
    tmp = tmp-0.5*id*gb*pow(0.5*(1-bout),(id-1));

  dmodeds = dmodeds+fa*tmp;

  /// Normalize
  dmodedr = pow(2,id+0.5)*dmodedr;
  dmodeds = pow(2,id+0.5)*dmodeds;
}

void rstoab(const fmatrix &r, const fmatrix &s,
	    fmatrix &a, fmatrix &b){

  ///function [a,b] = rstoab(r,s)
  /// Purpose : Transfer from (r,s) -> (a,b) coordinates in triangle

  int M = r.entryCount();
  a.resize(M,1);
  b.resize(M,1);
  for(int n=1;n<=M;++n){
    if(fabs(s(n)-1)>1e-5){
      a(n) = 2*(1+r(n))/(1-s(n))-1;
    }
    else
      a(n) = -1;
    b(n) = s(n);
  }
}

void Vandermonde2D(int p, fmatrix &rout, fmatrix &sout,
		   fmatrix &Vout){

  /// function [Vout] = Voutandermonde(N,r,s,Vout)
  /// Purpose : Initialize the gradient of the modal basis (i,j)
  ///	at (r,s) at order p

  Vout.resize(rout.entryCount(),(p+1)*(p+2)/2);

  /// find tensor-product coordinates
  fmatrix aout, bout;
  rstoab(rout,sout,aout,bout);

  /// Initialize matrices
  int sk = 1;
  for(int i=0;i<=p;++i){
    for(int j=0;j<=p-i;++j){
      for(int n=1;n<=rout.entryCount();++n){
	Vout(n,sk) = Simplex2DP(aout(n),bout(n),i,j);
      }
      sk = sk+1;
    }
  }
}


void GradVandermonde2D(int p, fmatrix &rout,
		       fmatrix &sout,
		       fmatrix &Vrout,
		       fmatrix &Vsout){

  /// function [Vrout,Vsout] = GradVandermonde2D(p,r,s)
  /// Purpose : Initialize the gradient of the modal basis (i,j)
  ///	at (r,s) at order N

  Vrout.resize(rout.entryCount(),(p+1)*(p+2)/2);
  Vsout.resize(rout.entryCount(),(p+1)*(p+2)/2);

  /// find tensor-product coordinates
  fmatrix aout, bout;
  rstoab(rout,sout,aout,bout);

  /// Initialize matrices
  int sk = 1;
  for(int i=0;i<=p;++i){
    for(int j=0;j<=p-i;++j){
      for(int n=1;n<=rout.entryCount();++n){
	GradSimplex2DP(aout(n),bout(n),i,j,Vrout(n,sk),Vsout(n,sk));
      }
      sk = sk+1;
    }
  }
}

void basis2D(int p, fmatrix &rout, fmatrix &sout,
	     fmatrix &Vout, fmatrix &Vrout,
	     fmatrix &Vsout){

  Vandermonde2D(p, rout, sout, Vout);
  GradVandermonde2D(p, rout, sout, Vrout, Vsout);
}

void Dmatrices2D(const int p, const fmatrix &r,
		 const fmatrix &s,
		 const fmatrix &V,
		 fmatrix &Dr,
		 fmatrix &Ds);

void Lift2D(fmatrix &Lift);
