SRCPATH=src
TCLAPINCLUDE=src/tclap-1.2.1/include
SOURCES=$(SRCPATH)/Util.cpp $(SRCPATH)/Bead.cpp $(SRCPATH)/Chromosome.cpp $(SRCPATH)/Randomizer.cpp $(SRCPATH)/Constraint.cpp $(SRCPATH)/Model.cpp $(SRCPATH)/MCMC.cpp $(SRCPATH)/Chrom3D.cpp

# Any C++11 compiler works (no boost needed), e.g. GCC 4.8.1 or newer, Clang 3.3 or newer, Apple clang.
# make uses its default C++ compiler (c++ or g++); choose another one with e.g.: make CXX=clang++   or   make CXX=g++-16
# Other options can be given in CXXFLAGS (default: -O3), e.g.: make CXXFLAGS="-O2 -Wall"
CXXFLAGS = -O3

# The C++ standard to compile as: the newest one the compiler supports, so that older compilers
# (which do not know e.g. -std=c++17) can build it too. The code is C++11, and every standard gives
# the same results. A standard is used only if the bundled TCLAP headers compile with it without
# warnings. Choose one yourself with e.g.: make CXXSTD=-std=c++11
ifndef CXXSTD
CXXSTD := $(shell for std in c++23 c++20 c++17 c++14 c++11; do \
            if echo '\#include <tclap/CmdLine.h>' | $(CXX) -std=$$std -Werror -I $(TCLAPINCLUDE) \
                 -x c++ -fsyntax-only - >/dev/null 2>&1; then echo -std=$$std; break; fi; done)
endif

# Always used:
# -ffp-contract=off: the compiler must not merge a*b+c into a fused multiply-add on its own. Compilers
#   do that differently (and it depends on the CPU and the optimisation level), which made the results
#   depend on the compiler. Where Chrom3D uses fused multiply-adds (on 64-bit ARM, e.g. Apple Silicon,
#   as the earlier versions did there) they are written out in the source (util::mulAdd, see
#   src/Util.h), so with this option all compilers give the same results on the same type of processor.
REQUIRED_FLAGS = $(CXXSTD) -ffp-contract=off

all: Chrom3D

Chrom3D: $(SOURCES) $(wildcard $(SRCPATH)/*.h) Makefile
	$(CXX) $(REQUIRED_FLAGS) $(CXXFLAGS) -I $(TCLAPINCLUDE) $(SOURCES) -o Chrom3D

clean:
	rm -f Chrom3D

.PHONY: all clean
