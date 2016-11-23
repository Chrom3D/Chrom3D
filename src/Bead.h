#ifndef GUARD_BEAD_H
#define GUARD_BEAD_H

#include <string>
#include <iostream>
#include <math.h>
#include "Util.h"


class Bead {
 public:
  Bead(int, int, double, std::string);
  Bead();
  void setCoordinates(double, double, double);  
  void setColor(float, float, float);
  double getX();
  double getY();
  double getZ();
  double getRadius();
  uint getBpSpan();
  std::string getID();
  bool hasCoordinates();
  
  friend std::ostream& operator<<(std::ostream&, const Bead&);
  friend bool operator<(const Bead&, const Bead&);
  friend double operator-(const Bead&, const Bead&);
  friend bool operator==(const Bead&, const Bead&);
  friend bool operator!=(const Bead&, const Bead&);      
  float color [3];
  void setRadius(double);
 private:
  double x;
  double y;
  double z;
  double radius;
  int start;
  int end;

  std::string id;
  bool hasCoord;
  int genomicPos;

  

};

bool beadClash(Bead&, Bead&);

#endif
