CFLAGS = -O3
SRCPATH=src
TCLAPINCLUDE=src/tclap-1.2.1/include

# Only a C++11 compiler is needed (no boost).
all:
	g++ -std=c++11 -O3 -I $(TCLAPINCLUDE) $(SRCPATH)/Util.cpp $(SRCPATH)/Bead.cpp $(SRCPATH)/Chromosome.cpp $(SRCPATH)/Randomizer.cpp $(SRCPATH)/Constraint.cpp $(SRCPATH)/Model.cpp $(SRCPATH)/MCMC.cpp $(SRCPATH)/Chrom3D.cpp -o Chrom3D

