#include "Randomizer.h"


using namespace std;

Move::Move(Chromosome &chr): chr(chr) {}

boost::container::static_vector<Bead, MAXSIZE>::iterator Move::begin() {
  return chr.begin() + start;
}

boost::container::static_vector<Bead, MAXSIZE>::iterator Move::end() {
  return chr.begin() + stop;
}

CoordinateMatrix::const_iterator1 Move::moveBegin() {
  return mat.begin1();
}

CoordinateMatrix::const_iterator1 Move::moveEnd() {
  return mat.end1();
}
 
void Move::validate() {
  assert(start <= stop);
  assert(mat.size1() == (stop-start)+1);
  assert(mat.size2() == 3);  
}

//////////////////////////////////////

Move CrankShaftRandomizer::randomize(Chromosome &chr, ENG &eng) {
  uint r1 = util::unif_int(0, chr.size()-3, eng);
  uint r2 = util::unif_int(r1+2, chr.size()-1, eng);
  double theta = util::unif_real(0.0, 2*util::PI-0.001, eng);
  CoordinateMatrix mat = chr.crankshaft(r1,r2,theta);
  Move mov(chr);
  mov.start = r1+1; // We return the positions of the start and end of the 
  mov.stop = r2-1;  //  beads that have their coords altered (so, r1 and r2 are not included in crankshaft)
  mov.mat = mat;
  return mov;
}

TranslateRandomizer::TranslateRandomizer(double mm){
  maxmove = mm;
}


Move TranslateRandomizer::randomize(Chromosome &chr, ENG &eng) {
   
  CoordinateVector tvec(3);
  tvec(0) = util::unif_int(-this->maxmove, this->maxmove, eng);
  tvec(1) = util::unif_int(-this->maxmove, this->maxmove, eng);
  tvec(2) = util::unif_int(-this->maxmove, this->maxmove, eng);

  CoordinateMatrix mat = chr.translate(tvec);
  Move mov(chr);
  mov.start = 0;
  mov.stop = chr.size()-1;
  mov.mat = mat;
  return mov;
}


Move ArmWiggleRandomizer::randomize(Chromosome &chr, ENG &eng) { 
  uint r1 = util::unif_int(1, chr.size()-2, eng);

  CoordinateVector origin(3);
  origin(0) = chr[r1].getX();
  origin(1) = chr[r1].getY();
  origin(2) = chr[r1].getZ();

  Move mov(chr);
  
  bool goRight = util::unif_int(0,1, eng);
  vector<double> radii;

  boost::container::static_vector<Bead,MAXSIZE>::iterator beaditer;
  if(goRight) {
    for(beaditer=chr.begin()+r1; beaditer!=chr.end(); beaditer++) {
      radii.push_back(beaditer->getRadius());
    }
    mov.start = r1;
    mov.stop = chr.size()-1;
      
  }
  else {
    for(beaditer=chr.begin(); beaditer!=chr.begin()+r1; beaditer++) {
      radii.push_back(beaditer->getRadius());
    }
    mov.start = 0;
    mov.stop = r1-1;            
  }

  CoordinateMatrix mat = util::createSelfAvoidingWalk(radii, origin, eng);
  mov.mat = mat;
  return mov;
}


Move ArmRotateRandomizer::randomize(Chromosome &chr, ENG &eng) {
  uint r1 = util::unif_int(1, chr.size()-2, eng);
  bool goRight = util::unif_int(0, 1, eng);
  double theta = util::unif_real(0.0, 2*util::PI-0.001, eng);
  
  CoordinateMatrix mat = chr.armRotate(r1, theta, goRight);
  Move mov(chr);
  if(goRight){
    mov.start = r1+1;
    mov.stop = chr.size()-1;
  } 
  else {
    mov.start = 0;
    mov.stop = r1-1;
  }
  mov.mat = mat;
  return mov;
}


Move MicroCrankShaftRandomizer::randomize(Chromosome &chr, ENG &eng) {
  uint r1 = util::unif_int(0, chr.size()-3, eng);
  uint r2 = r1+2;
  double theta = util::unif_real(0.0, 2*util::PI-0.001, eng);
  CoordinateMatrix mat = chr.crankshaft(r1,r2,theta);
  Move mov(chr);
  mov.start = r1+1; // We return the positions of the start and end of the 
  mov.stop = r2-1;  //  beads that have their coords altered (so, r1 and r2 are not included in crankshaft)
  mov.mat = mat;
  return mov;
}

Move RotationRandomizer::randomize(Chromosome &chr, ENG &eng) {
  uint r1 = util::unif_int(0, chr.size()-2, eng);
  uint r2 = util::unif_int(r1+1, chr.size()-1, eng);
  double theta = util::unif_real(0.0, 2*util::PI-0.001, eng);
  CoordinateMatrix mat = chr.rotation(r1,r2,theta);
  Move mov(chr);
  mov.start = 0; // We return the positions of the start and end of the 
  mov.stop = chr.size()-1;  //  beads that have their coords altered (so, r1 and r2 are not included in crankshaft)
  mov.mat = mat;
  return mov;


}
