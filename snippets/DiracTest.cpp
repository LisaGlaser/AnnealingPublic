#include "../TwoSphere.h"
#include "../progParams.h"
#include "../Dirac.h"
#include "../Randomc/stocc.h"
#include "../Randomc/randomc.h"
#include "../Randomc/rancombi.cpp"
#include "../Randomc/mersenne.cpp"
#include "../Randomc/mother.cpp"
#include "../Randomc/stoc1.cpp"                // random library source code
// define system specific user interface:
#include "../Randomc/userintf.cpp"

int main() {
	programParams iniV;
	iniV.initialconfig = 4;
	iniV.matrixsize = 30;
	iniV.Type = 00;

	Dirac D(iniV);
	std::vector <double> myevs = D.get_eigenvalues();
	std::cout << "eigenvalues: " << std::endl;
	for (auto i: myevs)
		std::cout << i << ' ';
	std::cout << std::endl;
	double ds = D.getds();
	std::cout << "spectral dimension: " << ds << std::endl;
}
