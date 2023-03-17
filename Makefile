CFLAGS = -O3
SRCPATH=src
TCLAPINCLUDE=src/tclap-1.2.1/include

all:
	g++ -Wall -Wextra -Wsign-compare -Wunused -Wpedantic -Wnull-dereference -O3 -I $(TCLAPINCLUDE) $(SRCPATH)/Util.cpp $(SRCPATH)/Bead.cpp $(SRCPATH)/Chromosome.cpp $(SRCPATH)/Randomizer.cpp $(SRCPATH)/Constraint.cpp $(SRCPATH)/Model.cpp $(SRCPATH)/MCMC.cpp $(SRCPATH)/Chrom3D.cpp -o Chrom3D

