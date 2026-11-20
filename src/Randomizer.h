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
  BeadIterator begin();
  BeadIterator end();
  void validate();
};


class Randomizer {
 public:
  virtual ~Randomizer() {} // The Model deletes its randomizers through Randomizer pointers
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
