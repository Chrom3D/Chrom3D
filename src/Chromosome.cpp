#include "Chromosome.h"
#include "Util.h"

#include <assert.h>
#include <stdexcept>


using namespace std;
namespace ublas = boost::numeric::ublas;
  
//using namespace boost::numeric::ublas;

Chromosome::Chromosome() {
  name = "";
  boost::container::static_vector<Bead,MAXSIZE> beadList;
  map<string,int> bead2index;
}


Chromosome::Chromosome(string chrName) {
  name = chrName;
  boost::container::static_vector<Bead,MAXSIZE> beadList;
  map<string,int> bead2index;
}

void Chromosome::addBead(Bead bead) {
  assert(bead.hasCoordinates());
  boost::container::static_vector<Bead,MAXSIZE>::iterator pos;
  pos = lower_bound(beadList.begin(), beadList.end(), bead);
  assert(pos == beadList.end() or bead < *pos); // Bead must not be found allready on the chr
  beadList.insert(pos, bead);  
}

void Chromosome::printBeads() {  
  for(boost::container::static_vector<Bead,MAXSIZE>::iterator it = beadList.begin(); it != beadList.end(); it++) {
    cout << *it << endl;
  }
}

uint Chromosome::size() {
  return this->beadList.size();
}

string Chromosome::getName() {
  return this->name;
}


boost::container::static_vector<Bead, MAXSIZE>::iterator Chromosome::begin() {
  return this->beadList.begin();
}

boost::container::static_vector<Bead, MAXSIZE>::iterator Chromosome::end() {
  return this->beadList.end();
}


CoordinateMatrix Chromosome::getBeadCoordinates(uint startPos, uint endPos) {
  assert(startPos <= endPos);
  assert(endPos < beadList.size());
  CoordinateMatrix mat((endPos-startPos)+1,3);
  uint i=0;
  for(boost::container::static_vector<Bead,MAXSIZE>::iterator it = beadList.begin() + startPos; it != beadList.begin() + endPos+1; it++) {
    mat(i,0) = it->getX();
    mat(i,1) = it->getY();
    mat(i,2) = it->getZ();
    i++;
  }
  return mat;
}

CoordinateVector Chromosome::getBeadCoordinates(uint pos) {
  return ublas::row(this->getBeadCoordinates(pos, pos),0);
}

CoordinateMatrix Chromosome::getBeadCoordinates() {
  return this->getBeadCoordinates(0, this->beadList.size()-1);
}


void Chromosome::printCoordinates(uint startPos, uint endPos) {
  cout << this->getBeadCoordinates(startPos, endPos);
}

void Chromosome::setColor(double r, double g, double b) {
  assert(r >= 0 and r <= 1);
  assert(g >= 0 and g <= 1);
  assert(b >= 0 and b <= 1);
  for(boost::container::static_vector<Bead,MAXSIZE>::iterator it = beadList.begin(); it != beadList.end(); it++) {
    it->setColor(r, g, b);
  }
}

CoordinateMatrix Chromosome::crankshaft(uint startPos, uint endPos, double theta) {
  // Add asserts here!
  CoordinateVector startPoint = this->getBeadCoordinates(startPos);
  CoordinateVector endPoint = this->getBeadCoordinates(endPos);
  CoordinateMatrix mat = this->getBeadCoordinates(startPos+1, endPos-1);
  CoordinateMatrix rotated = util::rotatePointsAroundArbitraryAxis(theta, startPoint, endPoint, mat);
  return rotated;
}

CoordinateMatrix Chromosome::rotation(uint startPos, uint endPos, double theta) {
  // Add asserts here!
  CoordinateVector startPoint = this->getBeadCoordinates(startPos);
  CoordinateVector endPoint = this->getBeadCoordinates(endPos);
  CoordinateMatrix mat = this->getBeadCoordinates(0, this->size()-1); // Full matrix is rotated
  CoordinateMatrix rotated = util::rotatePointsAroundArbitraryAxis(theta, startPoint, endPoint, mat);
  return rotated;
}



CoordinateMatrix Chromosome::armRotate(uint pos, double theta, bool goRight) {
  assert( (pos >= 1) && (pos < this->size()-1) );
  assert( (theta >= 0) && (theta <= 2*util::PI) );
  
  CoordinateVector beadPos = this->getBeadCoordinates(pos);
  CoordinateMatrix mat;
  CoordinateVector point1;
  CoordinateVector point2;
  if(goRight) {
    mat = this->getBeadCoordinates(pos+1, this->size()-1);
    point1 = this->getBeadCoordinates(pos-1);
    point2 = this->getBeadCoordinates(pos);
  }
  else {
    mat = this->getBeadCoordinates(0, pos-1);
    point1 = this->getBeadCoordinates(pos);
    point2 = this->getBeadCoordinates(pos+1);    
  }
  CoordinateMatrix rotated = util::rotatePointsAroundArbitraryAxis(theta, point1, point2, mat);  
  return rotated;
}

CoordinateMatrix Chromosome::translate(CoordinateVector translation) {
  return this->getBeadCoordinates() + util::repeatVector(translation,this->size());
}


bool operator==(const Chromosome& chr1, const Chromosome& chr2) {
  return chr1.name == chr2.name;
}

bool operator!=(const Chromosome& chr1, const Chromosome& chr2) {
  return !(chr1 == chr2);
}

Bead& Chromosome::operator[](const uint idx) {
  return beadList[idx];
}
