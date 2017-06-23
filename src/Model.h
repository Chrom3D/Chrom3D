#ifndef GUARD_MODEL_H
#define GUARD_MODEL_H

#include <vector>
#include "Chromosome.h"
#include "Bead.h"
#include "Constraint.h"
#include "Randomizer.h"
#include "Util.h"

#include <algorithm>
#include <string>
#include <map>
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/lexical_cast.hpp>
#include <iostream>
#include <fstream>
#include <numeric>

#include <boost/unordered_map.hpp>
#include <math.h>
#include <list>
#include <limits>

#include <boost/unordered_set.hpp>

#include <sstream>

typedef std::pair<std::string, std::string> idPair;

class Model {
 public:
  Model(std::string, uint seed=1234);
  std::string name;

  void addChromosome(Chromosome);  
  double getLossScore();
  double getLossScore(ConstraintType);
  double getPartialLossScore();  
  void addConstraint(std::string beadId1, std::string beadId2, double targetDistance, ConstraintType constType = REGULAR, double springConstant=1, double lowerBound=0, double upperBound=std::numeric_limits<double>::infinity());
  void addInteractionConstraint(std::string beadId1, std::string beadId2, double springConstant=1, ConstraintType interactionType = INTERACTION);
  void addInteractionDistanceConstraint(std::string beadId1, std::string beadId2, double beadPairDesiredDistance, double springConstant=1);
  void addInteractionLowerDistanceConstraint(std::string beadId1, std::string beadId2, double beadPairDesiredDistance, double springConstant=1);
  void addInteractionUpperDistanceConstraint(std::string beadId1, std::string beadId2, double beadPairDesiredDistance, double springConstant=1);
  void addBoundaryConstraint(std::string, double, double springConstant=1.0);
  void addNucleusConstraint(std::string, double springConstant=1.0);
  void addCenterConstraint(std::string, double springConstant=1.0);
  void addPeripheryConstraint(std::string, double springConstant=1.0);
  void addRandomizer(Randomizer*, double weight=1.0);
  void removeAllRandomizers();
  double getPartialLossScore(Chromosome&, uint, uint);
  Move proposeMove(uint maxAttempts=1);
  ENG randEngine;

  Move accept(Move& mov); 
  void writeCMM(std::string);
  std::string getCMM();
  void readGtrack(std::string, bool scaleBeadSizes=false, double nuclearOccupancy=0.2);
  void setNucleusRadius(double);
  double getTotalChrLength();
  uint numberOfChromosomes();
  void resetAllChromosomes(uint maxAttepts=1000);
  uint getNumberOfBeads();
  void addNucleusConstraints(double springConstant=1.0);
  void addCenterConstraints(double springConstant=1.0);
  void addBoundaryConstraints(double springConstant=1.0);
  void addSmartConstraints(double springConstant=1.0);

  
  void switchConstraints(ConstraintType, bool);
  void reweightConstraints(ConstraintType, double);
  std::vector<Chromosome>::iterator begin();
  std::vector<Chromosome>::iterator end();
  bool isClashing(Move&, bool onlyIntra=false);
  bool detectClash();

  bool isClashingFast(Move&);
  
  std::map<std::string, std::vector<Constraint> > constraints;
  std::map<std::string, std::vector<Constraint> > constraintsOFF;
  boost::ptr_vector<Randomizer> randomizers; //http://www.boost.org/doc/libs/1_55_0/libs/ptr_container/doc/tutorial.html
  std::vector<double> randomizerWeights;
  double getTotalBeadVolume();
  double getTotalGenomeLength();
  void rescaleBeadSizes(double);
  double getNuclearOccupancy();
  private: 

  std::vector<Chromosome> chromosomes; 
  bool hasNucleus;
  double nucleusRadius;

  std::map<std::string, Bead*> getBead;

  void updateBeadMap();
  void addConstraint(Constraint);
  std::vector<double> calculateScorePerConstraint(Chromosome&, uint, uint);
  Bead centerBead;
  bool interactionSpecifiedSymmetricly(std::vector<std::pair<std::string,std::string> > &interactionVector);
  bool interactionWeightsSpecifiedSymmetricly(std::map<idPair,double> &interactionWithWeight);

};
uint weightedSamplePos(std::vector<double>, ENG&);


#endif
