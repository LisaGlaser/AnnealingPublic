#ifndef __TWOSPHERE_H__
#define __TWOSPHERE_H__
// Class to implement properties of the two sphere. These are mainly needed to calculate the
// 1st order Heisenberg relation as the constraint in the main code.

// NB this class uses C++11 for lambda expressions. to port to C++98, just rewrite the lambda expressions as separate functions. C++11 is used by default in GCC 6.1 and later; in older versions, specify -std=c++11.

// this is copied from Dirac.h
//#include "Randomc/randomc.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cstring>
#include <math.h>
#include <Eigen/Sparse>
#include <Eigen/Eigenvalues>
#include <cmath>

// this is necessary for this class
#include <vector>
#include <algorithm>

// this is the type of the eigenbasis and individual eigenvectors |l,m>_{+-}
typedef std::vector<std::tuple<int,int,bool>> s2mode_vect;
typedef std::tuple<int,int,bool> s2mode;

class TwoSphere
{
public:
  // giving this an explicit constructor in hopes that it helps
  //TwoSphere();

  // sets the number of positive eigenvalues of D_{S^2} we will use
  void setsize(int);
  int getsize();
  // calculates coefficients from Dabrowski et al 2005; 2a = x+iy in Cartesian coords, 2b = 1-z.
  // the basis we use are the eigenspinors of D, see gracia-bondia et al ch. 9A.
  // l,m are doubled to be ints instead of halfints
  // that is, D|l,m>_\pm = \pm ((l/2) +1/2) |l,m>_\pm and m runs from l to -l in steps of 2
  // true means - and false means + so that \pm = (-1)^{sgn}.
  static double getSU2GeneratorACoefficients(int l, int m, bool sgn, int lp, int mp, bool sgnp);
  static double getSU2GeneratorBCoefficients(int l, int m, bool sgn, int lp, int mp, bool sgnp);
  // return the Bott projector p using the specified size.
  // the eigenvalues are sorted as follows: minus before plus (true before false), ascending in l, then ascending in m.
  Eigen::MatrixXcd getPMat();
  // return the Dirac of the given size
  Eigen::MatrixXcd getDiracS2();
  // return the pos part of the dirac
  Eigen::MatrixXcd getPosDiracS2();
  // return the sqrt of the pos part of the dirac
  Eigen::MatrixXcd getPosDiracS2sqrt();
  // double a matrix as A-> {{A,0}{0,A}}; this is applied to D to let L(H) act on H \otimes \C^2 and to take commutators of {{D,0},{0,D}} with PMat \in L(H \otimes \C^2)
  Eigen::MatrixXcd doubleMatrix(Eigen::MatrixXcd);
  // return a vector of all eigenspinors in the specified range
  s2mode_vect getEigenSpinors();
  // map {{A,B},{C,D}} to {A+D}
  Eigen::MatrixXcd RelativeTrace(Eigen::MatrixXcd);
  // return relativetrace((p-1/2)([D,p][D,p]))
  Eigen::MatrixXcd oneSidedEquationLHS(Eigen::MatrixXcd D, Eigen::MatrixXcd p);
  // return the gamma matrix on L(H)
  Eigen::MatrixXcd getGamma();
  // return ||<(p-1/2)[D,p][D,p]>-\gamma||_2^2
  double oneSidedEquationFrobSq(Eigen::MatrixXcd D, Eigen::MatrixXcd p);
private:
  // sizeS2 is the number of positive eigenvalues of the cutoff D
  int sizeS2;
  // spinorsize is dim(L^2(S^2,S)) times the number of positive eigenvalues of D: this is the dimension of the whole matrix D
  int spinorsize;
  // totalsize is dim (L^2(S^2,S) \otimes M_2(\C)), the dimension of the ambient space where the commutators are taken
  int totalsize;
  // this uses the function coeff of l,m,sgn,lp,mp,sgnp to construct a matrix
  Eigen::MatrixXcd setSpinorsizedMatrix(double (coeff)(int,int,bool,int,int,bool));
  // construct the matrices P,gamma,Diracs2 and the eigenspinors once
  Eigen::MatrixXcd produceLeftUpperPMat();
  Eigen::MatrixXcd produceRightUpperPMat();
  Eigen::MatrixXcd produceLeftLowerPMat();
  Eigen::MatrixXcd produceRightLowerPMat();
  Eigen::MatrixXcd producePMat();
  Eigen::MatrixXcd produceGamma();
  Eigen::MatrixXcd produceDiracS2();
  Eigen::MatrixXcd producePosDiracS2();

  s2mode_vect produceEigenSpinors();
  // cache PBott, gamma,eigenspinors
  Eigen::MatrixXcd PBott;
  Eigen::MatrixXcd gamma;
  Eigen::MatrixXcd DiracS2;
  Eigen::MatrixXcd PosDiracS2;
  s2mode_vect eigenSpinors;
  // the sorting function defines the ordering of the eigenspinor basis, for easy comparison with mathematica and nice-looking matrices
  static bool eigenSort(s2mode, s2mode);
};

#endif
