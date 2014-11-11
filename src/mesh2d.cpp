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
#include "mesh2d.hpp"

int compare_columns(const void *a, const void *b){
  int *aI = (int*) a;
  int *bI = (int*) b;

  int a1 = aI[2], a2 = aI[3];
  int b1 = bI[2], b2 = bI[3];

  if(b2>a2) return -1;
  if(a2>b2) return  1;

  if(b1>a1) return -1;
  if(a1>b1) return  1;

  return 0;

}

mesh2d::mesh2d(){
  Nverts = 0;
  Nfaces = 0;
  Nvertsperface = 0;
}

/// constructor: populates mesh with stuff from file named filename
mesh2d::mesh2d(setupAide &setup){

  int Nodes, allElements, Elements, dummy;

  Nverts = 3;
  Nfaces = 3;
  Nvertsperface = 2;

  datafloat minX = 1e9;
  datafloat maxX = -1e9;
  datafloat minY = 1e9;
  datafloat maxY = -1e9;

  string filename;

  if(!setup.getArgs("MESH NAME", filename)) throw -1;

  if(filename.find(".msh")!=string::npos){

    /// TW: need to read BCType from file
    ///     [ Wall, Outflow, Inflow, ... ]
    /// try - throw - catch for error handling
    try{
      /// open input file stream
      std::ifstream meshfile(filename.c_str());

      /// check to see if the meshfile is not open
      if(!meshfile.is_open()) throw -1;

      while(meshfile.good()){
	string line;
	getline(meshfile, line);

	if(line.find("$Nodes") != string::npos){

	  datafloat foo;

	  meshfile >> Nodes;
	  if(Nodes==0) throw -2;

	  /// create row vectors for the
	  /// x and y coordinates of the vertices
	  VX.resize(1, Nodes); /// exception handling inside matrix class
	  VY.resize(1, Nodes);

	  /// read vertex list 1:Nodes
	  for(int v=1;v<=Nodes;v=v+1){
	    /// scan line, ignoring vertex index
	    meshfile >> dummy >> VX(1,v) >> VY(1,v) >> dummy ;

	    minX = mymin(minX, VX(1,v));
	    maxX = mymax(maxX, VX(1,v));
	    minY = mymin(minY, VY(1,v));
	    maxY = mymax(maxY, VY(1,v));
	  }
	}

	if(line.find("$Elements") != string::npos){

	  /// read number of elements
	  meshfile >> allElements;
	  getline(meshfile, line);

	  /// record location in file
	  int marker = meshfile.tellg();
	  int etype, Tags;

	  Elements = 0;
	  for(int k=1;k<=allElements;k=k+1){

	    meshfile >> dummy >> etype >> Tags;
	    for(int t=1;t<=Tags;t=t+1)
	      meshfile >> dummy;
	    getline(meshfile, line); /// pull in rest of line from meshfile

	    if(etype==2) /// gmsh tri ?
	      Elements = Elements + 1;
	  }

	  /// rewind to beginning of element list
	  meshfile.seekg(marker);

	  /// shape EToV to hold Elements triples
	  EToV.resize(Elements, Nverts);
	  EToV = 0.;

	  /// scan in triangles
	  Elements = 1;
	  for(int k=1;k<=allElements;k=k+1){

	    meshfile >> dummy >> etype >> Tags;
	    for(int t=1;t<=Tags;t=t+1)
	      meshfile >> dummy;

	    if(etype==2){ /// gmsh tri ?
	      for(int v=1;v<=Nverts;++v)
		meshfile >> EToV(Elements,v) ;

	      Elements = Elements + 1;
	    }

	    /// oops needed this here in either case
	    getline(meshfile, line);
	  }
	  Elements = Elements-1;
	}
      }
    }
    catch(int err_code){
      if(err_code==-1) {
	std::cout << "ERROR: mesh(string filename) failed to open file" << std::endl;
	exit(-1); /// optionally exit
      }
    }
  }

  int newElements = 0;
  datafloat minJallowed = 1e-8;

  /// correct orientation of negative elements by permuting their vertices
  imatrix okflag(Elements,1);
  for(int k=1;k<=Elements;k=k+1){

    datafloat x1 = VX[EToV(k,1)], x2 = VX[EToV(k,2)], x3 = VX[EToV(k,3)];
    datafloat y1 = VY[EToV(k,1)], y2 = VY[EToV(k,2)], y3 = VY[EToV(k,3)];

    datafloat xr = 0.5*(x2-x1);
    datafloat xs = 0.5*(x3-x1);
    datafloat yr = 0.5*(y2-y1);
    datafloat ys = 0.5*(y3-y1);

    datafloat J = xr*ys-xs*yr;

    if(J<0){
      //// swap two vertices to force positive Jacobian
      int tmp = EToV(k,3);
      EToV(k,3) = EToV(k,2);
      EToV(k,2) = tmp;
      datafloat tmpx = x3, tmpy = y3;
      x3 = x2; y3 = y2;
      x2 = tmpx; y2 = tmpy;
    }

    xr = 0.5*(x2-x1);
    xs = 0.5*(x3-x1);
    yr = 0.5*(y2-y1);
    ys = 0.5*(y3-y1);

    J = xr*ys-xs*yr;

    datafloat lnx[3], lny[3];
    lnx[0] = yr,    lny[0] = -xr;
    lnx[1] = ys-yr, lny[1] = -xs+xr;
    lnx[2] =-ys,    lny[2] =  xs;
    int flag = 1;
    for(int f=1;f<=Nfaces;++f){
      datafloat lsJ = sqrt(lnx[f-1]*lnx[f-1]+lny[f-1]*lny[f-1]);
      datafloat lh  = J/lsJ;
      if(lh<1e-4)
	flag = 0;
    }
    okflag(k) = flag;
    if(flag){
      ++newElements;
    }
  }

  imatrix newEToV(newElements,3);

  //#define DEBUG_MESH
#ifdef DEBUG_MESH
  std::cout << " mesh2d::mesh2d reports newElements= " << newElements << endl;
#endif
  newElements = 0;
  datafloat Hmin = 1e9, Hmax = -1e9;
#if 1
  for(int k=1;k<=Elements;k=k+1){

    datafloat x1 = VX[EToV(k,1)], x2 = VX[EToV(k,2)], x3 = VX[EToV(k,3)];
    datafloat y1 = VY[EToV(k,1)], y2 = VY[EToV(k,2)], y3 = VY[EToV(k,3)];

    datafloat xr = 0.5*(x2-x1);
    datafloat xs = 0.5*(x3-x1);
    datafloat yr = 0.5*(y2-y1);
    datafloat ys = 0.5*(y3-y1);

    datafloat J = xr*ys-xs*yr;
    if(okflag(k)){

      ++newElements;

      newEToV(newElements,1) = EToV(k,1);
      newEToV(newElements,2) = EToV(k,2);
      newEToV(newElements,3) = EToV(k,3);

      datafloat lnx[3], lny[3];
      lnx[0] = yr,    lny[0] = -xr;
      lnx[1] = ys-yr, lny[1] = -xs+xr;
      lnx[2] =-ys,    lny[2] =  xs;
      int flag = 1;
      for(int f=1;f<=Nfaces;++f){
	datafloat lsJ = sqrt(lnx[f-1]*lnx[f-1]+lny[f-1]*lny[f-1]);
	datafloat lh  = J/lsJ;
	Hmax = mymax(Hmax, lh);
	Hmin = mymin(Hmin, lh);
      }

      //	cout << " done write " << endl;
    }
    if(J<0) std::cout << " warning: negative J " << std::endl;
  }
  Elements = newElements;
  EToV = newEToV;


#endif

#ifdef DEBUG_MESH
  std::cout << " newElements " << newElements << std::endl;
  std::cout << " max(mesh H) =  " << Hmax << std::endl;
  std::cout << " min(mesh H) =  " << Hmin << std::endl;
  std::cout << " mesh2d::mesh2d drops " << Elements - newElements << std::endl;
  std::cout << " size(EToV) = " << EToV.nrows() << ", " << EToV.ncolumns() << std::endl;
#endif

  std::cout << " mesh2d::mesh2d reports Elements= " << Elements << std::endl;
  std::cout << " minX = " << minX << std::endl;
  std::cout << " maxX = " << maxX << std::endl;
  std::cout << " minY = " << minY << std::endl;
  std::cout << " maxY = " << maxY << std::endl;


}

  /// find element to element connectivty
void mesh2d::connect(){

  int K = EToV.nrows();
  int Nv = VX.ncolumns();

  matrix <int> FToV(Nvertsperface+2, Nfaces*K);
  int vnums[3][2] = {{1,2}, {2,3}, {3,1}};
  FToV = 0.0;

  int counter = 1;
  for(int k=1;k<=K;k=k+1){
    for(int f=1;f<=Nfaces;f=f+1){
      int ns[Nvertsperface];
      FToV(1, counter) = k;
      FToV(2, counter) = f;

      for(int i=0;i<Nvertsperface;++i)
	ns[i] = EToV(k, vnums[f-1][i]);

      for(int i=0;i<Nvertsperface;++i){
	for(int j=i+1;j<Nvertsperface;++j){
	  if(ns[j]<ns[i]){
	    int tmp = ns[i];
	    ns[i] = ns[j];
	    ns[j] = tmp;
	  }
	}
      }

      for(int i=0;i<Nvertsperface;++i)
	FToV(3+i, counter) = ns[i];

      ++counter;
    }
  }

  /// sort by 3rd row (forgot column major convention)
  FToV.sort(compare_columns);

  /// populate EToE and EToF connectivity arrays
  /// 1-indexed
  EToE.resize(K, Nfaces); EToE = 0.0;
  EToF.resize(K, Nfaces); EToF = 0.0;
  BCType.resize(K, Nfaces); BCType = 0.0;

  /// now find neighbors
  for(counter=1;counter<=K*Nfaces-1;++counter){
    if(FToV(3,counter)==FToV(3,counter+1) &&
       FToV(4,counter)==FToV(4,counter+1)){

      int k1 = FToV(1,counter);
      int f1 = FToV(2,counter);
      int k2 = FToV(1,counter+1);
      int f2 = FToV(2,counter+1);

      /// matching
      EToE(k1,f1) = k2;
      EToE(k2,f2) = k1;
      EToF(k1,f1) = f2;
      EToF(k2,f2) = f1;
    }
  }

  for(int k=1;k<=K;k=k+1)
    for(int f=1;f<=Nfaces;f=f+1){
      if(EToE(k,f)==0){
	/// fix this later for correct boundary conditions
	BCType(k,f) = BC_Wall;
      }
    }
}


void mesh2d::diagnose(){
  datafloat x1, x2, x3, y1, y2, y3;
  datafloat x12, x13, x21, x23, x31, x32;
  datafloat y12, y13, y21, y23, y31, y32;
  datafloat l12, l13, l21, l23, l31, l32;

  matrix<datafloat > theta(3, 1);
  datafloat theta_max, theta_min;

  stringstream badElementCoordinates;

  int bad_elements = 0;
  const int K = EToV.nrows();
  /*! check for skewed elements */
  for(int k=1; k<=K; ++k){
    // extract the vertices
    x1 = VX(EToV(k, 1));      x2 = VX(EToV(k, 2));      x3 = VX(EToV(k, 3));
    y1 = VY(EToV(k, 1));      y2 = VY(EToV(k, 2));      y3 = VY(EToV(k, 3));

    x12 = x2 - x1;      x13 = x3 - x1;      x21 = -x12;
    x23 = x3 - x2;      x31 = -x13;      x32 = -x23;

    y12 = y2 - y1;      y13 = y3 - y1;      y21 = -y12;
    y23 = y3 - y2;      y31 = -y13;      y32 = -y23;

    l12 = sqrt(x12 * x12 + y12 * y12);
    l13 = sqrt(x13 * x13 + y13 * y13);
    l21 = l12;
    l23 = sqrt(x23 * x23 + y23 * y23);
    l31 = l13;
    l32 = l23;

    theta(1) = acos((x12 * x13 + y12 * y13) /(l12 * l13));
    theta(2) = acos((x21 * x23 + y21 * y23) /(l21 * l23));
    theta(3) = acos((x31 * x32 + y31 * y32) /(l31 * l32));

    if(theta.maxentry()/theta.minentry() > 6){
      std::cout<<"Warning : The ratio of max angle to min angle of element "<<k<<" is "<< theta.maxentry()/theta.minentry()<<std::endl;

      badElementCoordinates << " ( " << x1 << " , " << y1 << " )\n";

      bad_elements++;
    }

  }

  std::cout << "Total number of bad elements = " << bad_elements << std::endl;
  if(bad_elements > 0)
    std::cout << "Bad Elements:\n" << badElementCoordinates.str();
  //    assert(bad_elements == 0);
}

void mesh2d::debug(const char *name){

  /// debugging output for mesh
  std::ofstream os(name);

  os << "VX = "   << VX << std::endl;
  os << "VY = "   << VY << std::endl;
  os << "EToV = " << EToV << std::endl;
  os << "EToE = " << EToE << std::endl;
  os << "EToF = " << EToF << std::endl;

}
