/*
 * MCMC.cxx
 *
 * Copyright 2019 Lisa Glaser <glaser@univie.ac.at> and Abel Stern
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 *
 *
 */

// This is the main file of the code. It includes the generation of the Markov chain,
// as well as handling storing measurements.

#include "progParams.h"
#include "Dirac.h"
#include "Randomc/stocc.h"                     // define random library classes

#include <typeinfo>
#include <stdio.h>
#include <iostream>
#include <time.h> // guess what, this measures time!
#include <sys/time.h>
//#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Eigenvalues>

#include <highfive/H5File.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>

#undef H5_USE_BOOST
#define H5_USE_BOOST

#include <boost/multi_array.hpp>

#include "Randomc/randomc.h"
#include <csignal>

#ifndef MULTIFILE_PROJECT
// If compiled as a single file then include these cpp files,
// If compiled as a project then compile and link in these cpp files.
// Include code for the chosen random number generator:
#include "Randomc/rancombi.cpp"
#include "Randomc/mersenne.cpp"
#include "Randomc/mother.cpp"
#include "Randomc/stoc1.cpp"                // random library source code
// define system specific user interface:
#include "Randomc/userintf.cpp"
#endif

#define DEBUGMC 0
#define SUPERBUGMC 0
#define DEBUG	0
#define STEST 0
#define DDBUG 0


#define GAUSS 0 // uncommment this to get back to ordinary moves // this variable also exists in Dirac.h
//#define MOREEV 0

using Eigen::MatrixXcd;


using namespace Eigen;
using namespace std;

// I only define the random number generator once, so that successive calls give successive random numbers

TRandomCombined<CRandomMersenne,CRandomMother> RanGen(time(NULL));


sig_atomic_t stopFlag = 0;

void signalHandler( int signum) { // , Dirac D, programParams iniV) {
   cout << "Interrupt signal (" << signum << ") received.\n";
   stopFlag=1;

   // cleanup and close up stuff here
   // terminate program
   exit(signum);

}


// sets attributes of hdf5 data structures
void set_attributes(programParams iniV, HighFive::DataSet *set)
{
    /// now we add the attributes
    // what do we need?
    // constraints, matrix size, number of measurements type, cosmological constant
    string gitversion(string("hash=")+string(GIT_HASH)+string(", time=")+string(COMPILE_TIME)+string(", branch=")+string(GIT_BRANCH));
    HighFive::Attribute gitInfo=set->createAttribute<string>("git version info",
    HighFive::DataSpace::From(gitversion));
    gitInfo.write(gitversion);

    HighFive::Attribute ms = set->createAttribute<int>(
    "Matrix size", HighFive::DataSpace::From(iniV.matrixsize));
    ms.write(iniV.matrixsize);
    HighFive::Attribute sw = set->createAttribute<int>(
    "Number of sweeps", HighFive::DataSpace::From(iniV.stepnumber));
    sw.write(iniV.stepnumber);
    HighFive::Attribute ic = set->createAttribute<int>(
    "Initial configuration", HighFive::DataSpace::From(iniV.initialconfig));
    ic.write(iniV.initialconfig);
    HighFive::Attribute ty = set->createAttribute<int>(
    "Type", HighFive::DataSpace::From(iniV.Type));
    ty.write(iniV.Type);
    HighFive::Attribute vc = set->createAttribute<double>(
    "k factor", HighFive::DataSpace::From(iniV.kfac));
    vc.write(iniV.kfac);
    HighFive::Attribute Tini = set->createAttribute<double>(
    "initial Temperature", HighFive::DataSpace::From(iniV.T0));
    Tini.write(iniV.T0);
    HighFive::Attribute Tfin = set->createAttribute<double>(
    "final Temperature", HighFive::DataSpace::From(iniV.Tf));
    Tfin.write(iniV.Tf);
    HighFive::Attribute mr = set->createAttribute<double>(
    "Moves Radius", HighFive::DataSpace::From(iniV.wmoveA));
    mr.write(iniV.wmoveA);

}

// creates a hdf5 dataset for the Dirac operator
HighFive::DataSet create_setD(programParams iniV, char* name, int length, HighFive::File *file)
{
    // just a little check before I get started
    if(length<0 or iniV.stepnumber<0)
    {
        cout<<"Length or stepnumber are negative, you messed up!"<<endl;
        exit(1);
    }
    else
    {
    long unsigned int alength=length;
    long unsigned int steps=iniV.stepnumber;
    // first create the dataspace
    HighFive::DataSpace daspa = HighFive::DataSpace({steps,alength*alength});
    HighFive::DataSetCreateProps props;
    props.add(HighFive::Chunking(std::vector<hsize_t>{1, alength*alength}));

    // this should add compression
    int deflate_level=9;
    props.add(HighFive::Shuffle());
    props.add(HighFive::Deflate(deflate_level));

    HighFive::DataSet set = file->createDataSet(name, daspa, HighFive::AtomicType<complex<double>>(), props);

    set_attributes(iniV,&set);

    return set;
    }
}

// creates a hdf5 dataset
HighFive::DataSet create_set(programParams iniV, char* name, int length, HighFive::File *file)
{
 // just a little check before I get started
     if(length<0 or iniV.stepnumber<0)
    {
        cout<<"Length or stepnumber are negative, you messed up!"<<endl;
        exit(1);
    }
    else
    {
    long unsigned int steps=iniV.stepnumber;
    long unsigned int alength=length;
    // first create the dataspace
    HighFive::DataSpace daspa = HighFive::DataSpace({steps, alength}, {HighFive::DataSpace::UNLIMITED,alength});
    HighFive::DataSetCreateProps props;
    props.add(HighFive::Chunking(std::vector<hsize_t>{1, alength}));

    // this should add compression
    int deflate_level=9;
    props.add(HighFive::Shuffle());
    props.add(HighFive::Deflate(deflate_level));

    HighFive::DataSet set = file->createDataSet(name, daspa, HighFive::AtomicType<double>(), props);

    set_attributes(iniV,&set);

    return set;
}
}

// creates a hdf5 dataset with scalars
void scalar_create(HighFive::File *file, programParams iniV,string group_name)
{
    // create a dataset ready to contain vectors of the scalar observables
    vector <string> labels;
    //D.getR(),D.getcon(),D.getVol(), D.getds(),action(iniV,&D), 0, weightA
    labels.push_back("Curvature R");
    labels.push_back("CCM constraint");
    labels.push_back("Volume");
    labels.push_back("Spectral Dimension");
    labels.push_back("annealing Temperature");
    labels.push_back("Acceptance rate");

    int len=labels.size();

    vector<double> single_array;
    for(int i=0; i<len; i++)
    {
        HighFive::DataSpace daspa = HighFive::DataSpace({1,1}, {HighFive::DataSpace::UNLIMITED,1});
        HighFive::DataSetCreateProps props;
        props.add(HighFive::Chunking(std::vector<hsize_t>{1, 1}));
        // this should add compression
        int deflate_level=9;
        props.add(HighFive::Shuffle());
        props.add(HighFive::Deflate(deflate_level));

        HighFive::DataSet dataset = file->createDataSet(group_name+labels[i], daspa,HighFive::AtomicType<double>(),props);
        set_attributes(iniV,&dataset);
    }


}
// writes out the scalara dataset into the datafile
void scalar_writeout(HighFive::File *file, programParams iniV, vector<vector<double>> data, string group_name)
{
    // create a dataset ready to contain vectors of the scalar observables
    vector <string> labels;
    //D.getR(),D.getcon(),D.getVol(), D.getds(),action(iniV,&D), 0, weightA
    labels.push_back("Curvature R");
    labels.push_back("CCM constraint");
    labels.push_back("Volume");
    labels.push_back("Spectral Dimension");
    labels.push_back("annealing Temperature");
    labels.push_back("Acceptance rate");

    unsigned long int len=data[0].size();
    unsigned long int n_measure=data.size();

    vector<double> single_array;
    for(int i=0; i<len; i++)
    {
        single_array={};
        for(int j=0; j<n_measure;j++)
        {
            single_array.push_back(data[j][i]);
        }
        HighFive::DataSet dataset = file->getDataSet(group_name+labels[i]);
        dataset.resize({n_measure, 1});
        dataset.write(single_array);

    }


}

// main loop generating the MC loop
int MChain(programParams iniV)
{

    Dirac Dtemp(iniV);
    Dirac D(iniV);
    double S=D.getcon();
    double Stemp=S;
    // a bunch of variables I need for the thermal annealing algorithm
    double dCt=0,dSt=0,T=iniV.T0,Tfin=iniV.Tf,dC=0,T0=iniV.T0;
    long unsigned int truesizeLI=D.truesize;

    /// if my code uses an initial configuration I better get the right temperature for that, otherwise the won't work
    if(iniV.initialconfig ==2)
    {
      vector<double> temps;
      HighFive::File file(iniV.inifile, HighFive::File::ReadOnly);
      std::string DATASET_NAME("/best case/annealing Temperature");
      // we get the dataset
      HighFive::DataSet dataset = file.getDataSet(DATASET_NAME);
      dataset.read(temps);
      if (!temps.empty())
          T0=temps.back();
    }
    T=T0;

    signal(SIGINT, signalHandler);
    signal(SIGABRT, signalHandler);
    signal(SIGTERM, signalHandler);
    //	signal(SIGKILL, signalHandler); // can't catch this

    // Create a new file using the default property lists.

    HighFive::File file(iniV.outfile, HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate);


    vector <vector<double>> scalars;

    string scalar_group_name = "/scalars/";
    file.createGroup(scalar_group_name.c_str());

    scalar_create(&file,iniV,scalar_group_name);

    vector <vector<double>> best;

    string best_group_name = "/best case/";

    file.createGroup(best_group_name.c_str());
    scalar_create(&file,iniV,best_group_name);


    // create dataset for the eigenvalues
    HighFive::DataSet eigenvalues_set =create_set(iniV,"eigenvalues",truesizeLI,&file);
    // create dataset for the Dirac
    HighFive::DataSet dirac_set =create_setD(iniV,"Dirac",truesizeLI,&file);
    // save how the data is in there
    string docstring("This saves vector.push_back(D(x,y)) with x the outer and y the inner loop ");
    HighFive::Attribute saved = dirac_set.createAttribute<string>(
        "Data saving convention", HighFive::DataSpace::From(docstring));
    saved.write(docstring);

    cout<<"final matrix"<<endl;
    cout<<iniV.Type<<endl;
    string finalmat= "WTF?";
    if(iniV.Type==0)
    {
        finalmat = "H0";
    }
    if(iniV.Type==2)
    {
        finalmat = "H0, L0";

    }

    HighFive::DataSet P_dataset =
    file.createDataSet<complex<double>>(finalmat, HighFive::DataSpace::From(D.get_final()));
    set_attributes(iniV,&P_dataset);

    P_dataset.write(D.get_final());

    HighFive::DataSet Best_dataset =
    file.createDataSet<complex<double>>("best Dirac", HighFive::DataSpace::From(D.get_D()));
    set_attributes(iniV,&Best_dataset);

    Best_dataset.write(D.get_D());

    int safety=10; // flush to file every safety measurements, means we lose less on death
    int sweep,i=0,ttemp,pretemp=1;
    double Sbest=S;
    double weightA,weightM,p,ar=0.,ptemp,n1=0.,beta=1/T;


    // initialize best so that ti always contains something
    best.push_back({D.getR(),D.getcon(),D.getVol(), D.getds(),beta,0});

    Dtemp=D;

    sweep=10*iniV.matrixsize*iniV.matrixsize; // too much data, sweep size is increased to matrixsize**2

    weightA=iniV.wmoveA;


    int thermtime=100;
    // this function implements the monte carlo chain
    // first we 'burn in' and adjust the temperature
    while( (T>Tfin || i<thermtime*sweep)&&stopFlag==0)
    {
        i++;
        Dtemp=D;
        // then do a move

        Dtemp.moveA(weightA);


        // count moves1,
        n1+=1.;



        Stemp=Dtemp.getcon();
        dC=S-Stemp;
        ptemp=exp(beta*dC);
        p=RanGen.Random();
        if(DEBUGMC) printf("The old action is %g the new action is %g\n ",S,Stemp);

        if(SUPERBUGMC) printf("Move happens if %g > %g \n",ptemp, p);

        if(ptemp>p)
        {
            // for the tempering algorithm
            dCt-=dC;

            D=Dtemp;
            S=Stemp;
            ar+=1.;
            if(S<Sbest && i<thermtime*sweep)
            {
                Sbest=S;
                Best_dataset.write(D.get_D());
                best.push_back({D.getR(),D.getcon(),D.getVol(), D.getds(),T,ar/n1});
                //cout<<"New smallest CCM "<<Sbest<<endl;
            }


            if(DEBUGMC) printf("We got a move! \n");

        }
        if(dC>0)
        {
            dSt-=dC*beta;
        }

        // then measure the observables
        if(i%sweep== 0)
        {
            if(i%(sweep*safety)==0)
                {
                    // now the icing on the cake would be if I deleted scalars after every writeout, but I think the memory saving is not worth the trouble
                    P_dataset.write(D.get_final());
                    scalar_writeout(&file,iniV,best,best_group_name);
                    file.flush();
                }

            if(ttemp==1) // check if the measurement failed, if it did throw me out.
            {
                return 1;
            }

            double Tl=T;
            if(dCt>=0 || dSt==0) T=T0;
            else                 T=iniV.kfac*dCt/dSt;

            //weightA=weightA*sqrt(sqrt(T/Tl));
            //cout<<"new temperature "<<T<<" and new weight "<<weightA<<"Acceptance rate "<<ar<<" in step "<<i<<endl;
            //cout<<"dCt "<<dCt<<" and dSt "<<dSt<<endl;

            ar=0;
            n1=0;
            beta=1./T;

        }

    }
    printf("The final tempreature is %g which is below %g",T,iniV.Tf);
    // then we take the measurements
    i=0;

    while( i<iniV.stepnumber*sweep)
    {
        i++;
        Dtemp=D;
        // then do a move

        Dtemp.moveA(weightA);


        // count moves1,
        n1+=1.;

        Stemp=Dtemp.getcon();
        dC=S-Stemp;
        ptemp=exp(beta*dC);
        p=RanGen.Random();
        if(DEBUGMC) printf("The old action is %g the new action is %g\n ",S,Stemp);

        if(SUPERBUGMC) printf("Move happens if %g > %g \n",ptemp, p);

        if(ptemp>p)
        {

            D=Dtemp;
            S=Stemp;
            ar+=1.;
            if(S<Sbest)
            {
                Sbest=S;
                Best_dataset.write(D.get_D());
                best.push_back({D.getR(),D.getcon(),D.getVol(), D.getds(),T,ar/n1});
                //cout<<"New smallest CCM "<<Sbest<<endl;
            }

            if(DEBUGMC) printf("We got a move! \n");

        }

        // then measure the observables
        if(i%sweep== 0)
        {
            cout<<"measuring now "<<i/sweep<<"th time"<<endl;
            //	ttemp=measure(measureF,iniV,&D);
            long unsigned int row=i/sweep-1;
            // measuring the scalar things
            scalars.push_back({D.getR(),D.getcon(),D.getVol(), D.getds(),T,ar/n1});
            // eigenvalue writing
            eigenvalues_set.select({row,0},{1,truesizeLI}).write(D.get_eigenvalues());
            // dirac writing


            dirac_set.select({row,0},{1,truesizeLI*truesizeLI}).write(D.get_D());

            // getting the file to output data in between so I don't lose everything anymore.
            if(i%(sweep*safety)==0)
                {
                    // now the icing on the cake would be if I deleted scalars after every writeout, but I think the memory saving is not worth the trouble

                    P_dataset.write(D.get_final());
                    scalar_writeout(&file,iniV,scalars,scalar_group_name);
                    scalar_writeout(&file,iniV,best,best_group_name);
                    file.flush();
                }

            if(ttemp==1) // check if the measurement failed, if it did throw me out.
            {
                return 1;
            }


        }

    }

    Dtemp=D;
    scalar_writeout(&file,iniV,best,best_group_name);
    scalar_writeout(&file,iniV,scalars,scalar_group_name);

    file.flush();
    cout<<"end of chain"<<endl;
    return 0;
}

int main(int argc, char *argv[])
{
	// Some variables I will be using during the main loop


	programParams parameters;


    /// The first thing I do is initialize my program. For that I read a bunch of things from an input file

	// the programm expects the argument to be an input file
		if(argc==1)
		{
			cout<<"This programm needs an input file as an argument."<<endl;
			return 1;
		}
		else
		{

            /// A function to initialize everything
			parameters.initialize(argv[1]);

			parameters.announce();
			srand (time(NULL));

		cout<<"starting chain"<<endl;
			MChain(parameters);
		cout<<"end of code"<<endl;

        return 0;
      }
}
