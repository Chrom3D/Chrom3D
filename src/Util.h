#ifndef GUARD_UTIL_H
#define GUARD_UTIL_H

#include <stdint.h>
#include <numeric>
#include <vector>

//#include "Bead.h"
//#include "Chromosome.h"
#include <math.h> 
#include <iostream>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp> 
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/container/static_vector.hpp>



#include <boost/random.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <numeric>
#include <algorithm>

#include <map>

#include <exception>



typedef boost::numeric::ublas::matrix<double> CoordinateMatrix;
typedef boost::numeric::ublas::vector<double> CoordinateVector;
typedef unsigned int uint;

typedef boost::random::mt19937 ENG;

const double PRECISIONCONST=1.0e-15;

namespace ublas = boost::numeric::ublas;


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
  bool coordClash(CoordinateVector, CoordinateMatrix, double, std::vector<double>);
  CoordinateMatrix createSelfAvoidingWalk(std::vector<double>, CoordinateVector, ENG& eng, uint maxIter=1000);
  std::vector<std::string> &split(const std::string&, char, std::vector<std::string>&);
  std::vector<std::string> split(const std::string&, char);
  std::vector<double> &splitDbl(const std::string&, char, std::vector<double>&);
  std::vector<double> splitDbl(const std::string&, char);
  std::map<std::string, std::string> makeMap(std::vector<std::string>, std::vector<std::string>);

  class MoveException: public std::exception {
  public:
    MoveException();
    const char* what(){return "Unacceptable move!";}
  };

}
#endif
