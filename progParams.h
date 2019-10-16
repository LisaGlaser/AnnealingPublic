// This is start of the header guard.  initial_H can be any unique name.  By convention, we use the name of the header file.
#ifndef __PROGPARAMS_H__
#define __PROGPARAMS_H__

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cstring>
#include <math.h>

#define DEBUG 0
#define DBL_MAX  144115188075855872

// This class just defines a whole bunch of parameter values and the routines to read them in

class programParams
{

public:
	programParams();
	int initialize(char *filename);
	void announce();

	int matrixsize;
	int stepnumber;
	char *outfile;
	int initialconfig;
	char *inifile;
	double Lc;
	double T0;
	double Tf;
	double kfac;
	double wmoveA;
	int Type;

};


// This is the end of the header guard
#endif
