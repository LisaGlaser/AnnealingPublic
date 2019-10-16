/// A class that will carry all my parameters

#include "progParams.h"


programParams::programParams()
{

	// some hardcoded default values, for lazy people
	matrixsize=2;
	wmoveA=1./6./pow(matrixsize,1.5);
	stepnumber=1000;
	outfile= (char*) calloc (150,sizeof(char));
	strcpy(outfile,"default_output.txt");
	initialconfig=1;
	Type=00;
	Lc=-1;
	kfac=1;
	T0=100;
	Tf=0.01;
	inifile= (char*) calloc (150,sizeof(char));

}

int programParams::initialize(char *filename)
{
	FILE *fin;
	int j=0,iniI=0,buf;
	char * value;
	char buffer[150];
	double p;


	value= (char*) calloc (150,sizeof(char));
//	buffer= (char*) calloc (30,sizeof(char));

	fin = fopen (filename, "r");
	if (!fin)
	{
		printf("Sorry but your input file does not exist");
		return 1;
	}


	// while there are more lines we keep going
	while(j!=-1)
	{
		j=fscanf(fin, "%s  %s \n", value, buffer);

	if(strcmp(value,"matrixsize")==0)
	{

		buf=atoi(buffer);
		matrixsize=buf;
		wmoveA=1./6./pow(matrixsize,1.5);


		if(DEBUG) printf(" %d %s \n",buf,value);
	}
	else if(strcmp(value,"Type")==0)
	{	/// I am assuming that the type is (p,q) written as 9pq.
		/// this only works as long as p,q<10 but I think it's a safe bet
		buf=atoi(buffer);
		Type=buf;

		if(DEBUG) printf(" %d %s \n",buf, value);
	}
	else if(strcmp(value,"steps")==0)
	{
		buf=atoi(buffer);
		stepnumber=buf;

		if(DEBUG) printf(" %d %s \n",buf,value);
	}
	else if(strcmp(value,"outputfile")==0)
	{
		strcpy(outfile,buffer);

		if(DEBUG) printf(" %s %s \n",outfile,value);
	}
	else if(strcmp(value,"initialconfig")==0)
	{
		buf=atoi(buffer);
		initialconfig=buf;

		if(DEBUG) printf(" %d %s \n",buf,value);
	}
	else if(strcmp(value,"initalfile")==0)
	{
		strcpy(inifile,buffer);
		iniI=1;

		if(DEBUG) printf(" %s %s \n",inifile,value);
	}
	else if(strcmp(value,"movesradiusA")==0)
	{
		p=atof(buffer);
		wmoveA=p;

		 if(DEBUG) printf(" %e %s \n",p,value);
	}
	else if(strcmp(value,"cutoff")==0)
	{
		p=atof(buffer);
		Lc=p;

		 if(DEBUG) printf(" %e %s \n",p,value);
	}
	else if(strcmp(value,"factor")==0)
	{
		p=atof(buffer);
		kfac=p;

		 if(DEBUG) printf(" %e %s \n",p,value);
	}
	else if(strcmp(value,"Tini")==0)
	{
		p=atof(buffer);
		T0=p;

		 if(DEBUG) printf(" %e %s \n",p,value);
	}
	else if(strcmp(value,"Tfinal")==0)
	{
		p=atof(buffer);
		Tf=p;

		 if(DEBUG) printf(" %e %s \n",p,value);
	}
	else
	{
		printf("Are you sure %s is a valid option? \n",value);
	}

	}

	if(iniI==0&& initialconfig==4)
		{
			printf("You need to tell me an initial file.");
			return 1;
		}


	fclose(fin);
	free(value);
	//free(buffer);
	return 0;

}

void programParams::announce()
{

		std::cout<< "You are simulating a " <<matrixsize<< "x"<< matrixsize<< " matrix"<<std::endl;
		std::cout<< "Of type (p,q)=pq " <<Type<< " " <<std::endl;
		std::cout<< "There will be "<< stepnumber <<" sweeps of "<< matrixsize*4 << " attempted MC moves"<<std::endl;
		std::cout<< "The  cutoff  "<< Lc<<std::endl;
		std::cout<< "The  initial temperature is "<< T0<<std::endl;
		std::cout<< "and the final temperature is "<< Tf<<std::endl;
		if(wmoveA==0.){ std::cout<< " The additive move distance will be determined dynamically"<<std::endl; }
		else{		std::cout<< "And a additive move distance of "<< wmoveA <<std::endl;}

		std::cout<< "And the resulting data will be gathered in " << outfile <<std::endl;
}
