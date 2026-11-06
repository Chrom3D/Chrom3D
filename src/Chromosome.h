#ifndef GUARD_CHROMOSOME_H
#define GUARD_CHROMOSOME_H

#include <vector>
#include <map>
#include <string>
#include "Util.h"
#include "Bead.h"

typedef unsigned int uint;

// The beads of a chromosome are kept in a std::vector (before: a boost::container::static_vector
// with room for at most 5000 beads). Constraints keep pointers to beads, so beads must not be
// added to a chromosome after constraints are created (as with the static_vector).
typedef std::vector<Bead>::iterator BeadIterator;

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
  BeadIterator begin();
  BeadIterator end();

  Bead& operator[](const uint idx);

  friend bool operator==(const Chromosome&, const Chromosome&);
  friend bool operator!=(const Chromosome&, const Chromosome&);
  void setColor(double, double, double);
 private:
  bool fixed;
  std::string name;
  std::vector<Bead> beadList;

};



#endif
