#ifndef GUARD_CONSTRAINT_H
#define GUARD_CONSTRAINT_H


#include "Bead.h"
#include <limits>
#include <math.h>
#include <assert.h>
#include <utility> 

enum ConstraintType {REGULAR=0, INTERACTION=1, PERIPHERY=2, NUCLEUS=3, BOUNDARY=4, CENTER=5, C_T_COUNT=6, INTERACTION_INTRA=7, INTERACTION_INTER=8, INTERACTION_DIST=9, NON_INTERACTION_DIST=10};

class Constraint {
  public:
  Bead* bead1;
  Bead* bead2;
  ConstraintType constType;

  Constraint(Bead* b1, Bead* b2, double targetDistance, ConstraintType ct = REGULAR, double springConstant=1, double lowerBound=0, double upperBound=std::numeric_limits<double>::infinity());
  std::pair<std::string, std::string> getBeadIDs();
  double k;
  double lower;
  double upper;
  double distance;
  double eval();
  std::string getConstraintType();
  int getConstraintTypeID();
};

#endif


