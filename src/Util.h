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
    const char* what(){return "Unacceptable move!";}
  };

}
#endif
