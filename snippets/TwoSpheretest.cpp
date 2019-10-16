#include "../TwoSphere.h"
#include <iostream>


using namespace std;

int main()
{
  TwoSphere myTwoSphere;
  int size = 6;
  cout << "# Test for size " << size << endl;
  myTwoSphere.setsize(size);
  // cout<<myTwoSphere.getSU2GeneratorACoefficients(1,-1,1,1,1,1)<<endl;
  s2mode_vect eigenSpinors=myTwoSphere.getEigenSpinors();
  cout << "# sorted eigenvectors: eigenvalue, l, m, sign:" << endl;
  for (s2mode_vect::iterator i = eigenSpinors.begin(); i != eigenSpinors.end(); i++){
    int dl = get<0>(*i);
    int dm = get<1>(*i);
    int sgn = get<2>(*i);
    cout << (dl+1)/2 << " ";
    cout << dl << "/2 ";
    cout << dm << "/2 ";
    cout << pow(-1,sgn) << endl;
  }
  // cout << "This is the upper right block of P" << endl;
  // cout << myTwoSphere.getRightUpperPMat() << endl;
  cout << endl;
  cout << "# This is the P matrix" << endl;
  cout << myTwoSphere.getPMat() << endl;
  cout << endl;
  cout << "# This is the Dirac" << endl;
  cout << myTwoSphere.getDiracS2() << endl;
  cout << "# this is PosDiracS2" << endl;
  cout << myTwoSphere.getPosDiracS2() << endl;
  // int totalsize = 4*myTwoSphere.getsize();
  // Eigen::MatrixXcd D = myTwoSphere.DiracS2();
  // Eigen::MatrixXcd p = myTwoSphere.getPMat();
  // Eigen::MatrixXcd Dp = myTwoSphere.doubleMatrix(D);
  // cout << "This is double D" << endl << Dp << endl;
  // Eigen::MatrixXcd one = Eigen::MatrixXcd::Identity(totalsize,totalsize);
  // cout << "this is 1" << endl << one << endl;
  // Eigen::MatrixXcd commutator = Dp*p - p*Dp;
  // cout << "this is the comm" << endl << commutator << endl;
  // Eigen::MatrixXcd beforetrace=(p-one/2)*commutator*commutator;
  // cout << "this is beforetrace" << endl << beforetrace << endl;
  // cout << "This is he commutator identity" << endl;
  // cout << myTwoSphere.oneSidedEquationLHS(myTwoSphere.DiracS2(),myTwoSphere.getPMat()) << endl;
  cout << "# This is the LHS of the one-sided equation" << endl;
  cout << myTwoSphere.oneSidedEquationLHS(myTwoSphere.getDiracS2(),myTwoSphere.getPMat()) - myTwoSphere.getGamma() << endl;
  cout << "# This is the frob squared for the std dirac of size " << myTwoSphere.getsize() << endl;
  cout << myTwoSphere.oneSidedEquationFrobSq(myTwoSphere.getDiracS2(),myTwoSphere.getPMat()) << endl;
}

