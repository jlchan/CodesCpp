#include "midg1d.hpp"

static datafloat factorial(int n){

  if(n==0)
    return 1;
  else
    return n*factorial(n-1);

}


datafloat JacobiP(datafloat xout,
		  datafloat alpha,
		  datafloat beta, int p){

  // function [P] = JacobiP(x,alpha,beta,p)
  // Purpose: Evaluate Jacobi Polynomial of type (alpha,beta) > -1
  //          (alpha+beta <> -1) at points x for order N and
  //          returns P[1:length(xp))]
  // Note   : They are normalized to be orthonormal.

  // Turn points into row if needed.
  datafloat xp = xout;

  fmatrix PL(p+1,1);

  /// Initial values P_0(x) and P_1(x)
  datafloat gamma0 = pow(2,(alpha+beta+1))/(alpha+beta+1)*factorial(alpha)*factorial(beta)/factorial(alpha+beta);
  PL(1) = 1.0/sqrt(gamma0);
  if (p==0) return PL(1);

  datafloat gamma1 = (alpha+1)*(beta+1)/(alpha+beta+3)*gamma0;
  PL(2) = ((alpha+beta+2)*xp/2 + (alpha-beta)/2)/sqrt(gamma1);
  if (p==1) return PL(p+1);

  /// Repeat value in recurrence.
  datafloat aold = 2/(2+alpha+beta)*sqrt((alpha+1)*(beta+1)/(alpha+beta+3));

  /// Forward recurrence using the symmetry of the recurrence.
  for(int i=1;i<=p-1;++i){
    datafloat h1 = 2*i+alpha+beta;
    datafloat anew = 2/(h1+2)*sqrt( (i+1)*(i+1+alpha+beta)*(i+1+alpha)*(i+1+beta)/(h1+1)/(h1+3));
    datafloat bnew = -(alpha*alpha-beta*beta)/h1/(h1+2);
    PL(i+2) = 1./anew*( -aold*PL(i) + (xp-bnew)*PL(i+1));
    aold =anew;
  }

  return PL(p+1);
}

datafloat GradJacobiP(datafloat xout, datafloat alpha,
		      datafloat beta,  int p){

  /// function [dP] = GradJacobiP(z, alpha, beta, N);
  /// Purpose: Evaluate the derivative of the orthonormal Jacobi
  ///	   polynomial of type (alpha,beta)>-1, at points x
  ///          for order N and returns dP[1:length(xp))]

  datafloat dP;
  if(p == 0)
    dP = 0.0;
  else
    dP = sqrt(p*(p+alpha+beta+1))*JacobiP(xout,alpha+1,beta+1, p-1);

  return dP;
}
