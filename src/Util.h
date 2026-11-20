#ifndef GUARD_UTIL_H
#define GUARD_UTIL_H

#include <stdint.h>
#include <numeric>
#include <vector>
#include <string>
#include <sstream>

//#include "Bead.h"
//#include "Chromosome.h"
#include <math.h>
#include <cmath>
#include <iostream>
#include <random>
#include <numeric>
#include <algorithm>

#include <map>

#include <exception>
#include <stdexcept>

#include "Coordinates.h"


typedef unsigned int uint;

typedef std::mt19937 ENG; // Same generator (and same numbers for a given seed) as boost::random::mt19937

const double PRECISIONCONST=1.0e-15;


namespace util {
  const double PI=3.14159265358979323846264338327950;

  // x*x. Used instead of pow(x,2): when optimising, Clang replaces pow(x,2) by x*x (so that is what
  // Chrom3D computed), but the pow function itself sometimes differs from x*x in the last bit, and not
  // every compiler or optimisation level makes that replacement.
  inline double square(double x) { return x * x; }

  // a*b+c, rounded as the earlier versions of Chrom3D (built with their Makefile) rounded it on this
  // type of processor: on 64-bit ARM (e.g. Apple Silicon) as one fused multiply-add (a single rounding),
  // because Clang fuses a*b+c written in one expression there; on other processors (e.g. x86-64, whose
  // basic instruction set has no fused multiply-add) as a multiplication followed by an addition.
  // The build turns off any other fusing by the compiler (-ffp-contract=off, see the Makefile), so the
  // results depend on the type of processor (as before) but not on the compiler or optimisation level.
#if defined(__aarch64__) || defined(__arm64__)
  inline double mulAdd(double a, double b, double c) { return std::fma(a, b, c); }
#else
  inline double mulAdd(double a, double b, double c) { return a * b + c; }
#endif

  // Sets s = sin(theta) and c = cos(theta) in rotatePointsAroundArbitraryAxis. On macOS, Clang (when
  // optimising) computed the two together there with the math library's __sincos_stret, which sometimes
  // differs in the last bit from separate sin() and cos() calls. It is called directly so that the result
  // does not depend on the compiler. (createSelfAvoidingWalk uses separate calls, see Util.cpp.) Other
  // systems use sin() and cos(), which compilers merge where they did before (GCC into sincos).
  inline void sinCos(double theta, double &s, double &c) {
#if defined(__APPLE__)
    struct __double2 res = __sincos_stret(theta);
    s = res.__sinval;
    c = res.__cosval;
#else
    s = sin(theta);
    c = cos(theta);
#endif
  }

  double euclideanDistance(CoordinateVector, CoordinateVector);
  CoordinateMatrix repeatVector(CoordinateVector, uint);
  CoordinateVector normalizeVector(CoordinateVector);
  CoordinateMatrix rotatePointsAroundArbitraryAxis(double, CoordinateVector, CoordinateVector, CoordinateMatrix);
  double unif_real(double, double, ENG&);
  int unif_int(uint, uint, ENG&);
  bool coordClash(double, double, double, double, double, double, double, double);
  bool coordClash(CoordinateVector, CoordinateVector, double, double);
  bool coordClash(CoordinateVector, const CoordinateMatrix&, uint, double, const std::vector<double>&);
  CoordinateMatrix createSelfAvoidingWalk(std::vector<double>, CoordinateVector, ENG& eng, uint maxIter=1000);
  std::vector<std::string> &split(const std::string&, char, std::vector<std::string>&);
  std::vector<std::string> split(const std::string&, char);
  std::vector<double> &splitDbl(const std::string&, char, std::vector<double>&);
  std::vector<double> splitDbl(const std::string&, char);
  std::map<std::string, std::string> makeMap(std::vector<std::string>, std::vector<std::string>);
  std::string errorLine(uint line_counter, std::string filename);

  // Converts a whole string to a number, e.g. parseNumber<double>("0.2"). Like boost::lexical_cast
  // (which it replaces), a string that is not entirely a number, such as "1,0.5", is an error.
  template<typename T> T parseNumber(const std::string &s) {
    std::istringstream iss(s);
    T value;
    if(!(iss >> value) or !(iss >> std::ws).eof()) {
      throw std::invalid_argument("Could not read '" + s + "' as a number");
    }
    return value;
  }

  class MoveException: public std::exception {
  public:
    MoveException();
    const char* what() const noexcept override {return "Unacceptable move!";}
  };

}
#endif
