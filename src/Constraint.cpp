#include "Constraint.h"
using namespace std;

Constraint::Constraint(Bead* b1, Bead* b2, double targetDistance, ConstraintType ct/* = REGULAR */, double springConstant /*=1*/, double lowerBound /*=0*/, double upperBound /*=INF*/) {
  bead1 = b1;
  bead2 = b2;
  k = springConstant;
  distance = targetDistance;
  constType = ct;
  lower = lowerBound;
  upper = upperBound;
}


pair<string, string> Constraint::getBeadIDs() {
  return make_pair(bead1->getID(), bead2->getID());
}


double Constraint::eval() {
  double beadDist=(*bead1) - (*bead2);
  if(beadDist > this->lower and beadDist < this->upper) {
    return this->k * pow(beadDist - this->distance,2);
  }
  return 0;
}

string Constraint::getConstraintType(){

  switch(this->constType){
    case INTERACTION: return "INTERACTION";
    case PERIPHERY: return "PERIPHERY";
    case NUCLEUS: return "NUCLEUS";
    case BOUNDARY: return "BOUNDARY";
    case CENTER: return "CENTER";
    case INTERACTION_INTRA: return "INTERACTION_INTRA";
    case INTERACTION_INTER: return "INTERACTION_INTER";
    case INTERACTION_DIST: return "INTERACTION_DIST";
    case REGULAR: return "REGULAR";
    case C_T_COUNT: return "C_T_COUNT";
   }
   return "REGULAR";
}

int Constraint::getConstraintTypeID(){
  int id = this->constType;
  return id;
}
//DistanceConstraint::DistanceConstraint(Bead* b1, Bead* b2, double springConstant /*=1*/) : Constraint(b1,b2,b1->getRadius()+b2->getRadius(), springConstant) {};

