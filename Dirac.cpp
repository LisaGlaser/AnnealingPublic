#include "Dirac.h"


Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es;

Dirac::Dirac(programParams iniV)
{
	type=iniV.Type;
	size=iniV.matrixsize;
	if(type==2)
		{
			truesize=2*size;
		}
	else if(type==0)
		{
			truesize=2*size;
		}

  	// the dimension of D is twice the argument supplied here, so 2*size
	S2.setsize(size);

	initial(iniV.initialconfig,iniV.inifile);
	makeD();
	es.compute(D,false); // computes the eigenvalues, false tells it not to compute the eigenvectors
	eivals=es.eigenvalues();
	LcI=iniV.Lc; // I eventually want to remove Lc from the code completely
	Lc=LcI;
	LcFactors();
	setgamma5();
	calc_RVd();
	calc_S2con();
	moved=0;

}

Dirac::~Dirac()
{
	char finale[250];
	std::cout<<"Deconstructing Dirac"<<std::endl;
	strcpy(finale,"deconstructor_finalmatrix.txt");
	//printAll(finale);

	// This desctructor should also clean up my entire memory
	/// and I need to catch the sigkill signal to make it work
	//delete D;

}

void Dirac::dynamicLc()
{
	// no need to initialize, any square is larger than a negative number
	// unless we get imaginary numbers here, but then we are screwed anyway
	for(int i=0;i<eivals.size();i++)
	{
		if(eivals[i]*eivals[i]>Lc) Lc=eivals[i]*eivals[i];
	}
 //number from Mathematica
}

void Dirac::LcFactors()
{
	if(LcI<0) dynamicLc();
	E10=log(Lc+1)/(Lc+.001);
	// numbers fixed on 06/06/18 hopefully good now
	fLd=-75.39822368615503772310344120; // number copied from Mathematica sign ok?
	volumecoeff=21.0722*E10/2; //Number from mathematica
	// factor of 2 is the trace of the clifford module
}


double Dirac::F20(double l)
{ // we have fixed the dimension to 2
  return exp(-l*E10)/(1. + l*E10);
}

double Dirac::F21(double l)
{

	return exp(-4.*l*E10)/(1. +4.*l*E10) - exp(- 2.*l*E10)/(1. + 2.*l*E10);
}

void Dirac::calc_RVd()
{
	/// calculates R and V both because they go together so well
	double tunedM = 0.39612512715394627; // m is tuned so that ds = 2.0 for the sphere at size 30
	double S=0,rpart=0,volpart=0,kt=0,ktl2=0,ex=0;
	double ep=tunedM*log(Lc+1)*E10; // where 0.15 is the abel special factor

	// We got the normal eigenvalues and can just square them
	// and then throw away eigenvalues larger than Lc

	if(LcI<0) LcFactors(); /// actually not so good, now it does not recompute

	for(int x=0;x<eivals.size();x++)
	{
		ex=eivals[x];

		rpart+=F21(ex*ex);
	    volpart+=F20(ex*ex);

		kt+=exp(-1.*ex*ex*ep);
		ktl2+=ex*ex*exp(-1.*ex*ex*ep);

	}

	ds=2*ep*ktl2/kt; // get the spectral dimension by evaluating t p(t,D^2) / p(t,D) at t=\epsilon
	R=fLd*rpart;
	Vol=volumecoeff*volpart;
}

void Dirac::calc_S2con()
{
  //assume D is of dimension spinorsize, so already of the form {{-D_+,0},{0,D_+}}
  //we use S2.getPMat() now (which just returns the private Eigen::MatrixXcd S2.PBott calculated upon setsize() in the program init), because P is not yet dynamic.
 con = S2.oneSidedEquationFrobSq(D, S2.getPMat());

}


void Dirac::initial(int ini,char *inifile)
{

	Eigen::MatrixXcd mtemp(size,size);

	if(ini==1) // an identity Dirac operator
	{
		mtemp=Eigen::MatrixXcd::Identity(size, size);

		H0=mtemp;
		L0=Eigen::MatrixXcd::Zero(size,size);

	}
	else if(ini==2) // initial config from a file
	{
		H0.resize(size,size);
		L0.resize(size,size);

		std::vector<std::complex<double>> vecH0;
		HighFive::File file(inifile, HighFive::File::ReadOnly);
		std::string DATASET_NAME("WTF?");
		if(type==0)
		{
			DATASET_NAME="H0";
			L0=Eigen::MatrixXcd::Zero(size,size);
		}
		else if(type==2)
		{
			DATASET_NAME="H0, L0";
		}

		// we get the dataset
		HighFive::DataSet dataset = file.getDataSet(DATASET_NAME);

		// we convert the hdf5 dataset to a single dimension vector
		dataset.read(vecH0);


		read_H0(vecH0);
		std::cout<<H0<<std::endl;
		std::cout<<L0<<std::endl;


	}
	else if(ini==3) // random Dirac operator
	{

		mtemp=Eigen::MatrixXcd::Random(size,size);
		H0=mtemp+mtemp.adjoint();
		mtemp=Eigen::MatrixXcd::Random(size,size);
		L0=mtemp-mtemp.adjoint();

	}
  	else if(ini==4) // positive part of the S2 dirac
  	{
      H0 = S2.getPosDiracS2();
	  if(type==0)
	  	{
		H0 = S2.getPosDiracS2sqrt();
		//std::cout<<"type 0 get a sqrt!"<<std::endl;
		}
	  L0=Eigen::MatrixXcd::Zero(size,size);
	 // std::cout<<"initialized to"<<std::endl;
      //std::cout<<H0<<std::endl;
  	}
	else
	{
		std::cout<< " Sorry we don't support that many initial configurations yet"<<std::endl;
	}


}

std::vector<std::complex<double>> Dirac::get_D()
{

	std::vector<std::complex<double>> Dvec;

	for(int x=0;x<truesize;x++)
	{
		for(int y=0;y<truesize;y++)
		{
			Dvec.push_back(D(x,y));
		}
	}
	return Dvec;
}

std::vector<std::complex<double>> Dirac::get_final()
{

	std::vector<std::complex<double>> Dvec;

	if(type==0)
	{
		for(int x=0;x<size;x++)
		{
			for(int y=0;y<size;y++)
			{
				Dvec.push_back(H0(x,y));
			}
		}

	}
	else if(type==2)
	{
		for(int x=0;x<size;x++)
		{
			for(int y=0;y<size;y++)
			{
				Dvec.push_back(L0(x,y));
			}
		}
		for(int x=0;x<size;x++)
		{
			for(int y=0;y<size;y++)
			{
				Dvec.push_back(H0(x,y));
			}
		}

	}

	return Dvec;
}

void Dirac::read_H0(std::vector<std::complex<double>> H0v)
{


	if(type==0)
	{
		for(int x=0;x<size;x++)
		{
			for(int y=0;y<size;y++)
			{
				//Dvec.push_back(H0(x,y));
				H0(x,y)=H0v[x*size+y];
			}
		}

	}
	else if(type==2)
	{
		for(int x=0;x<size;x++)
		{
			for(int y=0;y<size;y++)
			{
				//Dvec.push_back(L0(x,y));
				L0(x,y)=H0v[x*size+y];
			}
		}
		for(int x=0;x<size;x++)
		{
			for(int y=0;y<size;y++)
			{
				H0(x,y)=H0v[size*size+x*size+y];
			}
		}

	}

	//return Dvec;
}

std::vector <double> Dirac::get_eigenvalues()
{
	std::vector <double> evs;

	for(int x=0;x<eivals.size();x++)
	{
	 	evs.push_back(eivals[x]);

	 }

	return evs;
}

double Dirac::getR()
{
	if(moved==1)
	{
		if(type==00) // only the upper corner is dynamical, hence we only need the eigenvalues for half the spectrum the other half is the negatives of it
			{es.compute(P,false); // computes the eigenvalues, false tells it not to compute the eigenvectors
			eivals<< es.eigenvalues(),-es.eigenvalues(); }
		else
			{
			es.compute(D,false); // computes the eigenvalues, false tells it not to compute the eigenvectors
			eivals=es.eigenvalues(); }
		calc_RVd();
		calc_S2con();
	}
	moved=0;
	return R;
}

double Dirac::getVol()
{
	if(moved==1)
	{
		if(type==00) // only the upper corner is dynamical, hence we only need the eigenvalues for half the spectrum the other half is the negatives of it
			{es.compute(P,false); // computes the eigenvalues, false tells it not to compute the eigenvectors
			eivals<< es.eigenvalues(),-es.eigenvalues(); }
		else
			{
			es.compute(D,false); // computes the eigenvalues, false tells it not to compute the eigenvectors
			eivals=es.eigenvalues(); }
		calc_RVd();
		calc_S2con();
	}
	moved=0;
	return Vol;
}


double Dirac::getds()
{
	if(moved==1)
	{
		if(type==00) // only the upper corner is dynamical, hence we only need the eigenvalues for half the spectrum the other half is the negatives of it
			{es.compute(P,false); // computes the eigenvalues, false tells it not to compute the eigenvectors
			eivals<< es.eigenvalues(),-es.eigenvalues(); }
		else
			{
			es.compute(D,false); // computes the eigenvalues, false tells it not to compute the eigenvectors
			eivals=es.eigenvalues(); }
		calc_RVd();
		calc_S2con();
	}
	moved=0;
	return ds;
}


double Dirac::getcon()
{
	if(moved==1)
	{
		if(type==00) // only the upper corner is dynamical, hence we only need the eigenvalues for half the spectrum the other half is the negatives of it
			{es.compute(P,false); // computes the eigenvalues, false tells it not to compute the eigenvectors
			eivals<< es.eigenvalues(),-es.eigenvalues(); }
		else
			{
			es.compute(D,false); // computes the eigenvalues, false tells it not to compute the eigenvectors
			eivals=es.eigenvalues(); }
		calc_RVd();
		calc_S2con();
	}
	moved=0;
	return con;
}

/// obviously right now this is hardcoded but if we have a better idea what we want we can make gamma5 dynamic
void Dirac::setgamma5()
{
	//std::cout<<size<<std::endl;
	gamma5=Eigen::MatrixXcd::Zero(size, size);

	for(int x=0;x<size;x++) gamma5.real()(x,size-x-1)=1;
}


void Dirac::moveA(double pA)
{
	int hermitian;
	Eigen::MatrixXcd temp;

	hermitian=1;

	moveA_raw(H0,hermitian,pA,size);
	if(!H0.isApprox(H0.adjoint())) std::cout<<"H0 warning"<<std::endl;
	temp=H0+H0.adjoint();
	H0=0.5*temp;

	hermitian=-1;

	moveA_raw(L0,hermitian,pA,size);
	if(!L0.isApprox(-L0.adjoint())) std::cout<<"L0 warning"<<std::endl;
	temp=L0-L0.adjoint();
	L0=0.5*temp;

	makeD();
	moved=1;
	if(!D.isApprox(D.adjoint()))
	{
		std::cout<<"That Dirac is not selfadjoint after moveM, something is wrong here!"<<std::endl;
		std::cout<<D-D.adjoint()<<std::endl;
	}

}



int seed = (int)time(0);    // random seed
StochasticLib1 sto(seed);           // make instance of random library



int moveA_raw(Eigen::MatrixXcd &M, int her, double p,int s)
{
	Eigen::MatrixXcd diff(s,s);
	diff=Eigen::MatrixXcd::Random(s,s);

	M+=sto.Normal(p,p)*(diff+ her*diff.adjoint()); // here comes the anti hermitian


	return 0;

}

void Dirac::makeD()
{


	if(type==00) // this is our default starting point, fixing the topology
	{
		makeD00();
	}
	else if(type==2) // this will be implemented later, no fixed topology, hence we will need to use two matrices, one hermitian one anti hermitian
	{
		makeDm2();
	}
	else
	{
		printf("Sorry, we haven't implemented Type %d",type);
	}


}

// Topology fixed, only a positive hermitian matrix.
// this does kind of disturb our measure though
void Dirac::makeD00()
{
	// we are only calculating the positive corner, and keep the dirac that way so truesize= size
	// need to finish fix trueszie first
	D=Eigen::MatrixXcd::Zero(truesize,truesize);

	P = H0*H0.adjoint();
	D.topLeftCorner(size,size) = -P;
	/// gamma5 is the permutation operator

	D.bottomRightCorner(size,size)=P;


}
// completely free, one hermitian and one anti hermitian matrix
void Dirac::makeDm2()
{
	// need to finish fix trueszie first
	D=Eigen::MatrixXcd::Zero(truesize,truesize);



	D.topLeftCorner(size,size) = H0;
	D.topRightCorner(size,size) = L0;
	/// gamma5 is the permutation operator

	D.bottomRightCorner(size,size)=-H0;
	D.bottomLeftCorner(size,size)=-L0;

}
