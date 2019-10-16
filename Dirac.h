#ifndef __DIRAC_H__
#define __DIRAC_H__

#include "progParams.h"
//#include "Randomc/randomc.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cstring>
#include <math.h>
#include <Eigen/Sparse>
#include <Eigen/Eigenvalues>
#include <cmath>
#include <highfive/H5File.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>

#undef H5_USE_BOOST
#define H5_USE_BOOST

#include <boost/multi_array.hpp>

#include "Randomc/stocc.h"                     // define random library classes
#include "TwoSphere.h"
//#define CUTOFF 0


// This class defines the Dirac operator, this includes how the moves act on the Dirac operator
// and how several observables are calculated from the Dirac operator

int moveA_raw(Eigen::MatrixXcd &M, int her, double p, int s);
int moveT_raw(Eigen::MatrixXcd &M, double p, int s);


class Dirac
{
	public:
//		Dirac();
		Dirac(programParams iniV);
		~Dirac();
		int size;
		int truesize;
		int type;

		double getR();
		double getVol();
		double getcon();
		double getds();

		std::vector <double> get_eigenvalues(); // turns my eigenvalues into a vector
		std::vector<std::complex<double>> get_D();// turns my Dirac into a complex array
		std::vector<std::complex<double>> get_final();// turns my Dirac into a complex array
		void moveA(double pA);
		void moveT(double pA);
		void makeD();

	private:
		double Lc,fLd,confac,volfac,volumecoeff,LcI; /// High energy cutoff lambda, a number in the action, a factor arising in the constraint, a factor in the volume, is the cutoff lambda dynamic

		double R; // curvature of the geometry
		double Vol; // volume of the geomtry
		double con; // constraint
		double ds; // spectral dimension

		double E10; // number in the action
		double hardcutoff; // no eigenvalue is allowed above Lambda

		Eigen::VectorXd eivals;
		Eigen::MatrixXcd D;
		Eigen::MatrixXcd H0,P,L0; // P is the positive matrix obtained as H0*H0.adjoint(), H0 is hermitian, L0 is antihermitian
		Eigen::MatrixXcd gamma5;

		TwoSphere S2;

		int moved;

		void read_H0(std::vector<std::complex<double>> H0v);
		void setU();
		void setgamma5();

		void calc_S2con();
		void calc_RVd();
		void dynamicLc();
		void LcFactors();
		double F21(double l);
    	double F20(double l);

		void eigenvalues();

		void makeD00(); // topology fixed
		void makeDm2(); // free

		void initial(int i,char *inifile);


};

#endif
