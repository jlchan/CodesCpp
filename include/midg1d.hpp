#ifndef __MIDG1D
#define __MIDG1D

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
#include "matrix.hpp"

datafloat JacobiP(datafloat xout, datafloat alpha,
		  datafloat beta, int p);

datafloat GradJacobiP(datafloat xout, datafloat alpha,
		      datafloat beta,  int p);

void JacobiGL(const datafloat alpha, const datafloat beta,
	      const int N, fmatrix &x);

void JacobiGQ(const datafloat alpha, const datafloat beta,
	      const int N, fmatrix &x, fmatrix &w);


#endif
