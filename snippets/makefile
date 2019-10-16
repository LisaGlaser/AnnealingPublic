# I am a comment, and I want to say that the variable CC will be
# the compiler to use.
CC=g++
#CC=clang++
# I'm using CFlags for compiler options, and just assume that LDFLAGS works for libraries to include. I might be abusing it here
CFLAGS= -c -g3 -std=c++11  -I../Eigen/Eigen3.24 -I../highfive/include -I/usr/include/hdf5/serial -Wno-ignored-attributes
#LDFLAGS= -lhdf5_serial ### needed on ubuntu
LDFLAGS= -lhdf5 ### should work otherwise

GIT_HASH=`git rev-parse HEAD`
COMPILE_TIME=`date -u +'%Y-%m-%d %H:%M:%S UTC'`
GIT_BRANCH=`git branch | grep "^\*" | sed 's/^..//'`
export VERSION_FLAGS=-DGIT_HASH="\"$(GIT_HASH)\"" -DCOMPILE_TIME="\"$(COMPILE_TIME)\"" -DGIT_BRANCH="\"$(GIT_BRANCH)\""


all: snippets

snippets: progParams.o Dirac.o TwoSphere.o DiracTest.o TwoSpheretest.o
	$(CC) progParams.o Dirac.o TwoSphere.o DiracTest.o -o DiracTest $(LDFLAGS)
	$(CC) TwoSphere.o TwoSpheretest.o -o TwoSpheretest $(LDFLAGS)

MCMC.o: ../MCMC.cpp
	$(CC) $(CFLAGS) ../MCMC.cpp $(LDFLAGS) $(VERSION_FLAGS)

progParams.o: ../progParams.cpp
	$(CC) $(CFLAGS) ../progParams.cpp

Dirac.o: ../Dirac.cpp
	$(CC) $(CFLAGS)  ../Dirac.cpp

TwoSphere.o: ../TwoSphere.cpp
	$(CC) $(CFLAGS) ../TwoSphere.cpp

DiracTest.o: DiracTest.cpp
	$(CC) $(CFLAGS) DiracTest.cpp

TwoSpheretest.o: TwoSpheretest.cpp
	$(CC) $(CFLAGS) TwoSpheretest.cpp

clean:
	rm *.o DiracTest
