#include "Bead.h"

#include <assert.h>

using namespace std;
Bead::Bead() {
  start = 0;
  end = 0;
  radius  = 0;
  id = "NA";
  genomicPos = 0;
  hasCoord = false;
  color[0] = 0;
  color[1] = 0;
  color[2] = 0;  
}

Bead::Bead(int startPosition, int endPosition, double beadRadius, std::string beadId) {
  assert(startPosition < endPosition);
  
  start = startPosition;
  end = endPosition;  
  radius = beadRadius;
  id = beadId;

  genomicPos = int((start + end)/2);
  hasCoord = false;

  color[0] = 1;
  color[1] = 1;
  color[2] = 0;
  
}

void Bead::setCoordinates(double xCoord, double yCoord, double zCoord) {
  x = xCoord;
  y = yCoord;
  z = zCoord;
  hasCoord = true;  
}

void Bead::setRadius(double newRadius) {
  this->radius = newRadius;
}

double Bead::getX() {
  return this->x;
}

double Bead::getY() {
  return this->y;
}

double Bead::getZ() {
  return this->z;
}


double Bead::getRadius() {
  return this->radius;
}

uint Bead::getBpSpan() {
  return(this->end - this->start);
}

bool Bead::hasCoordinates() {
  return hasCoord;
}

string Bead::getID() {
  return this->id;
}


void Bead::setColor(float red, float green, float blue) {
  color[0] = red;
  color[1] = green;
  color[2] = blue;
}


ostream& operator<<(ostream& os, const Bead& b) {
  os << b.id;
  return os;
}

double operator-(const Bead& b1, const Bead& b2) {
  return sqrt(pow(b1.x-b2.x,2) + pow(b1.y-b2.y,2) + pow(b1.z-b2.z,2));
}

bool operator<(const Bead& b1, const Bead& b2) {
  return b1.genomicPos < b2.genomicPos;
}

bool operator==(const Bead& b1, const Bead& b2) {
  return b1.id == b2.id;
}

bool operator!=(const Bead& b1, const Bead& b2) {
  return !(b1 == b2);
}
