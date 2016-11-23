#ifndef GUARD_CHROMOSOME_H
#define GUARD_CHROMOSOME_H

#include <vector>
#include <map>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp> 
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/container/static_vector.hpp>
#include <string>
#include "Util.h"
#include "Bead.h"

//typedef boost::container::static_vector<Bead,1000>::iterator beadIter;
typedef boost::numeric::ublas::matrix<double> CoordinateMatrix;
typedef boost::numeric::ublas::vector<double> CoordinateVector;
typedef unsigned int uint;

const uint MAXSIZE=1000; // Maximum size (number of beads) for each chromosome

class Chromosome {
 public:
  Chromosome(std::string);
  Chromosome();
  void addBead(Bead);
  void printBeads();
  void printCoordinates(uint, uint);
  std::string getName();
  uint size();
  CoordinateMatrix crankshaft(uint, uint, double);
  CoordinateMatrix rotation(uint, uint, double);
  CoordinateMatrix armRotate(uint, double, bool);
  CoordinateMatrix translate(CoordinateVector);
  CoordinateMatrix getBeadCoordinates(uint, uint);
  CoordinateVector getBeadCoordinates(uint);
  CoordinateMatrix getBeadCoordinates();
  boost::container::static_vector<Bead, MAXSIZE>::iterator begin();
  boost::container::static_vector<Bead, MAXSIZE>::iterator end();

  Bead& operator[](const uint idx);

  friend bool operator==(const Chromosome&, const Chromosome&);
  friend bool operator!=(const Chromosome&, const Chromosome&);
  void setColor(double, double, double);
 private:
  bool fixed;
  std::string name;
  //std::vector<Bead> beadList;
  boost::container::static_vector<Bead,MAXSIZE> beadList;
  
};



#endif
