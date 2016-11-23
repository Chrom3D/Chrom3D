#ifndef GUARD_RANDOMIZER_H
#define GUARD_RANDOMIZER_H

#include "Chromosome.h"
#include "Util.h"



class Move {
 public:
  Move(Chromosome&);
  Chromosome &chr;
  uint start;
  uint stop;
  CoordinateMatrix mat;
  boost::container::static_vector<Bead, MAXSIZE>::iterator begin();
  boost::container::static_vector<Bead, MAXSIZE>::iterator end();
  CoordinateMatrix::const_iterator1 moveBegin();
  CoordinateMatrix::const_iterator1 moveEnd();
  void validate();
};


class Randomizer {
 public:
  virtual Move randomize(Chromosome&, ENG&) = 0;
};


class CrankShaftRandomizer: public Randomizer {
 public:
  Move randomize(Chromosome&, ENG&);
};

class MicroCrankShaftRandomizer: public Randomizer {
 public:
  Move randomize(Chromosome&, ENG&);
};

class TranslateRandomizer: public Randomizer {
 public:
  TranslateRandomizer(double);
  Move randomize(Chromosome&, ENG&);
  
 private:
  double maxmove;
};

class ArmWiggleRandomizer: public Randomizer {
 public:
  Move randomize(Chromosome&, ENG&);
};

class ArmRotateRandomizer: public Randomizer {
 public:
  Move randomize(Chromosome&, ENG&);
};

class RotationRandomizer: public Randomizer {
 public:
  Move randomize(Chromosome&, ENG&);
};

#endif
