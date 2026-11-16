SRCPATH=src
TCLAPINCLUDE=src/tclap-1.2.1/include
SOURCES=$(SRCPATH)/Util.cpp $(SRCPATH)/Bead.cpp $(SRCPATH)/Chromosome.cpp $(SRCPATH)/Randomizer.cpp $(SRCPATH)/Constraint.cpp $(SRCPATH)/Model.cpp $(SRCPATH)/MCMC.cpp $(SRCPATH)/Chrom3D.cpp

# Any C++17 compiler works (no boost needed), e.g. GCC 8 or newer, Clang 7 or newer, Apple clang.
# make uses its default C++ compiler (c++ or g++); choose another one with e.g.: make CXX=clang++   or   make CXX=g++-16
# Other options can be given in CXXFLAGS (default: -O3), e.g.: make CXXFLAGS="-O2 -Wall"
CXXFLAGS = -O3

# Always used:
# -ffp-contract=off: the compiler must not merge a*b+c into a fused multiply-add on its own. Compilers
#   do that differently (and it depends on the CPU and the optimisation level), which made the results
#   depend on the compiler. Where Chrom3D uses fused multiply-adds they are written out in the source
#   (std::fma, see src/Util.cpp), so with this option all compilers give the same results.
REQUIRED_FLAGS = -std=c++17 -ffp-contract=off

all: Chrom3D

Chrom3D: $(SOURCES) $(wildcard $(SRCPATH)/*.h) Makefile
	$(CXX) $(REQUIRED_FLAGS) $(CXXFLAGS) -I $(TCLAPINCLUDE) $(SOURCES) -o Chrom3D

clean:
	rm -f Chrom3D

.PHONY: all clean
