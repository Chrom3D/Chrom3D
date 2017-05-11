#include "Model.h"
#include <assert.h>
#include <typeinfo>


using namespace std;

Model::Model(string modelName, uint seed /*=1234*/) {
  name=modelName;
  randEngine.seed(seed);
  centerBead.setCoordinates(0,0,0);
  hasNucleus=false;

}
vector<Chromosome>::iterator Model::begin() {
  return chromosomes.begin();
}

vector<Chromosome>::iterator Model::end() {
  return chromosomes.end();
}



void Model::updateBeadMap() {
  this->getBead.clear();
  vector<Chromosome>::iterator chriter;
  boost::container::static_vector<Bead,MAXSIZE>::iterator beaditer; 
  for(chriter = chromosomes.begin(); chriter != chromosomes.end(); chriter++) {
    for(beaditer = chriter->begin(); beaditer != chriter->end(); beaditer++) {
      getBead[beaditer->getID()] = &(*beaditer);
    }
  }
}


void Model::addChromosome(Chromosome chr) {
  assert ( find(chromosomes.begin(), chromosomes.end(), chr) == chromosomes.end() ); // Chromosomes must be unique
  // Assert that all Bead IDS ARE UNIQUE HERE AND NOT FOUND ON ANY CHR HERE
  chromosomes.push_back(chr);  
  this->updateBeadMap();
  //hasUpdatedHashMap = false;
}

void Model::addConstraint(Constraint constraint) {
  pair<string, string> beadIds = constraint.getBeadIDs();
  constraints[beadIds.first].push_back(constraint);
  //constraints[beadIds.second].push_back(constraint);  
}

void Model::addConstraint(string beadId1, string beadId2, double targetDistance, ConstraintType ct /*=REGULAR*/, double springConstant /*=1*/, double lowerBound /*=0*/, double upperBound /*=INF*/) {
  assert(getBead.find(beadId1) != getBead.end()); 
  assert(getBead.find(beadId2) != getBead.end()); 
  Bead* b1 = getBead[beadId1];
  Bead* b2 = getBead[beadId2];
  this->addConstraint(Constraint(b1, b2, targetDistance, ct, springConstant, lowerBound, upperBound));
}

void Model::addInteractionConstraint(string beadId1, string beadId2, double springConstant /*=1*/, ConstraintType interactionType/*= INTERACTION*/) {
  assert(getBead.find(beadId1) != getBead.end()); 
  assert(getBead.find(beadId2) != getBead.end());
  assert(interactionType == INTERACTION || interactionType == INTERACTION_INTRA || interactionType == INTERACTION_INTER);

  Bead* b1 = getBead[beadId1];
  Bead* b2 = getBead[beadId2];
  this->addConstraint(beadId1, beadId2, b1->getRadius()+b2->getRadius(), interactionType, springConstant);
}

void Model::addInteractionDistanceConstraint(string beadId1, string beadId2, double beadPairDesiredDistance, double springConstant /*=1*/, ConstraintType interactionType/*= INTERACTION_DIST*/) {
  assert(getBead.find(beadId1) != getBead.end()); 
  assert(getBead.find(beadId2) != getBead.end());
  assert(interactionType == INTERACTION_DIST);

  Bead* b1 = getBead[beadId1];
  Bead* b2 = getBead[beadId2];
  this->addConstraint(beadId1, beadId2, b1->getRadius()+b2->getRadius() + beadPairDesiredDistance, interactionType, springConstant, b1->getRadius()+b2->getRadius() + beadPairDesiredDistance);
}


void Model::addNonInteractionDistanceConstraint(string beadId1, string beadId2, double beadPairRequiredAboveDistance, double springConstant /*=1*/, ConstraintType interactionType/*= INTERACTION_DIST*/) {
  assert(getBead.find(beadId1) != getBead.end()); 
  assert(getBead.find(beadId2) != getBead.end());
  assert(interactionType == NON_INTERACTION_DIST);

  Bead* b1 = getBead[beadId1];
  Bead* b2 = getBead[beadId2];
  this->addConstraint(beadId1, beadId2, b1->getRadius()+b2->getRadius() + beadPairRequiredAboveDistance, interactionType, springConstant, 0, b1->getRadius()+b2->getRadius() + beadPairRequiredAboveDistance);
}

// Method addBoundary Constraint is not used at the moment
void Model::addBoundaryConstraint(string beadId, double boundaryRadius, double springConstant /*=1*/) {
  assert(getBead.find(beadId) != getBead.end()); 
  Bead* b1 = getBead[beadId];  
  this->addConstraint(Constraint(b1, &(this->centerBead), boundaryRadius - b1->getRadius(), BOUNDARY, springConstant, boundaryRadius - b1->getRadius()));
}


void Model::addNucleusConstraint(string beadId, double springConstant /*=1*/) {
  assert(getBead.find(beadId) != getBead.end()); 
  assert(this->hasNucleus);
  Bead* b1 = getBead[beadId];  
  this->addConstraint(Constraint(b1, &(this->centerBead), this->nucleusRadius - b1->getRadius(), NUCLEUS, springConstant, this->nucleusRadius - b1->getRadius()));
}

void Model::addLaminConstraint(string beadId, double springConstant /*=1*/) {
  assert(this->hasNucleus);
  assert(getBead.find(beadId) != getBead.end()); 
  Bead* b1 = getBead[beadId];  
  this->addConstraint(Constraint(b1, &(this->centerBead), this->nucleusRadius - b1->getRadius(), LAMIN, springConstant));
}

void Model::addCenterConstraint(string beadId, double springConstant /*=1*/) {
  assert(getBead.find(beadId) != getBead.end()); 
  Bead* b1 = getBead[beadId];  
  this->addConstraint(Constraint(b1, &(this->centerBead), 0, CENTER, springConstant));
}

void Model::setNucleusRadius(double r) {
  assert(r > 0);
  nucleusRadius = r;
  hasNucleus = true;
}

void Model::reweightConstraints(ConstraintType ct, double newSpringConstant) {
  map<string,vector<Constraint> >::iterator cm_it;
  vector<Constraint>::iterator c_it;   
  for(cm_it = constraints.begin(); cm_it != constraints.end(); cm_it++) {
      for(c_it = cm_it->second.begin(); c_it != cm_it->second.end();c_it++) {
        if(c_it->constType == ct) {
	  c_it->k = newSpringConstant;
      }
    }
  }
}

void Model::switchConstraints(ConstraintType ct, bool switchOn) {
  
  map<string,vector<Constraint> >::iterator cm_it;
  vector<Constraint>::iterator c_it; 
  pair<string, string> beadIds;
  
  if(switchOn){

    for(cm_it = constraintsOFF.begin(); cm_it != constraintsOFF.end(); cm_it++) {

      for(c_it = cm_it->second.begin(); c_it != cm_it->second.end();) {

        if(c_it->constType == ct) {
          beadIds = c_it->getBeadIDs();
          // First, copy all constraint pointers to the "ON" map (NOTE! We should check whether the constraint is already in the list) 
          constraints[beadIds.first].push_back(*c_it);
          // Erase the constraint from the "OFF" map vector for that bead
          cm_it->second.erase(c_it);
        }
        else {
          c_it++;
        }
      }
    }
  }
  else {

    for(cm_it = constraints.begin(); cm_it != constraints.end(); cm_it++) {

      for(c_it = cm_it->second.begin(); c_it != cm_it->second.end();) {

        if(c_it->constType == ct) {
          beadIds = c_it->getBeadIDs();
          // First, copy all constraint pointers to the "OFF" map (NOTE! We should check whether the constraint is already in the list) 
          constraintsOFF[beadIds.first].push_back(*c_it);
          // Erase the constraint from the "OFF" map vector for that bead
          cm_it->second.erase(c_it);
        }
        else {
          c_it++;
        }
      }
    }
  }
}


double Model::getLossScore() {
  map<string,vector<Constraint> >::iterator it1;
  vector<Constraint>::iterator it2;  
  double res = 0;
  string constraintType;
  for(it1 = constraints.begin(); it1 != constraints.end(); it1++) {
      for(it2 = it1->second.begin(); it2 != it1->second.end(); it2++) {
	constraintType = it2->getConstraintType();
	if(constraintType == "INTERACTION" || constraintType == "INTERACTION_DIST" || constraintType == "NON_INTERACTION_DIST" || constraintType == "INTERACTION_INTRA" || "INTERACTION_INTER") {
	  res += it2->eval()/2;
	}
	else {
	  res += it2->eval();
	}
      }
  }
  return res; 
}

double Model::getLossScore(ConstraintType ctype) {
  map<string,vector<Constraint> >::iterator it1;
  vector<Constraint>::iterator it2;  
  double res = 0;
  string constraintType;
  for(it1 = constraints.begin(); it1 != constraints.end(); it1++) {
      for(it2 = it1->second.begin(); it2 != it1->second.end(); it2++) {
	if(it2->constType == ctype) {
	  constraintType = it2->getConstraintType();
	  if(constraintType == "INTERACTION" || constraintType == "INTERACTION_DIST" || constraintType == "NON_INTERACTION_DIST" || constraintType == "INTERACTION_INTRA" || "INTERACTION_INTER") {
	   res += it2->eval()/2;
	  }
	  else {
	    res += it2->eval();
	  }
	}
      }
  }
  return res; 
}




double Model::getPartialLossScore(Chromosome &chr, uint startPos, uint endPos) {
  boost::container::static_vector<Bead,MAXSIZE>::iterator it;
  vector<Constraint>::iterator it2;
  double res = 0;
  for(it = chr.begin() + startPos; it != chr.begin() + endPos+1; it++) {
    for(it2 = constraints[it->getID()].begin(); it2 != constraints[it->getID()].end(); it2++) {
      // Within a local evaluation of the loss score, any pair of beads with interactions between them (inside the move), will be counted twice. Therefore, we need to divide (by 2) the loss-score of each constraint belonging to such  pairs. (Note that e.g. Nucleus constraints will never be down-weighted like this, since the centerBead is never inside the local move, and that is the intended behavior, since there is only one of those):
      if ((find(chr.begin()+startPos, chr.begin()+endPos+1, *(it2->bead1)) != chr.begin()+endPos+1) and (find(chr.begin()+startPos, chr.begin()+endPos+1, *(it2->bead2)) != chr.begin()+endPos+1) ) {
	res += it2->eval() / 2;
      }
      else {
	res += it2->eval();
      }
    }
  }
  return res;
}

void Model::addRandomizer(Randomizer* randomizer, double weight/*=1.0*/) {
  randomizers.push_back(randomizer);
  randomizerWeights.push_back(weight);
}

void Model::removeAllRandomizers() {
  
  randomizers.clear();
  
  randomizerWeights.clear();
}


Move Model::proposeMove(uint maxAttempts/*=10000*/) {
  assert(this->randomizers.size() > 0);
  uint iRand = weightedSamplePos(randomizerWeights, randEngine);
  uint iChr = util::unif_int(0, chromosomes.size()-1, randEngine);
  for(uint i=0; i!=maxAttempts; i++) {
    Move mov = randomizers[iRand].randomize(chromosomes[iChr], randEngine);
    if(mov.mat.size1()!=0) {
      mov.validate();
      if (not isClashing(mov)) {
        return mov;
      }
    }
  }
  throw util::MoveException();

}


bool Model::isClashing(Move &mov, bool onlyIntra /*=false*/) { //

  vector<Chromosome>::iterator chriter;  
  boost::container::static_vector<Bead,MAXSIZE>::iterator beaditer;

  uint movSize = (mov.stop-mov.start)+1;
  for(uint i=0; i!=movSize; i++) {
    double moveX=mov.mat(i,0);
    double moveY=mov.mat(i,1);
    double moveZ=mov.mat(i,2);
    double moveRadius=mov.chr[mov.start+i].getRadius();

    // INTRACHROMOSOMAL
    // "left" part:
    for(beaditer = mov.chr.begin(); beaditer != mov.begin(); beaditer++) {
      if(util::coordClash(beaditer->getX(),beaditer->getY(), beaditer->getZ(), moveX, moveY, moveZ, beaditer->getRadius(), moveRadius)) {
        return true;
      }
    }
    // "Right" part:
    for(beaditer = mov.end()+1; beaditer != mov.chr.end(); beaditer++) {
      if(util::coordClash(beaditer->getX(),beaditer->getY(), beaditer->getZ(), moveX, moveY, moveZ, beaditer->getRadius(), moveRadius)) {
        return true;
      }
    }
  }
    
  if(onlyIntra) {
    return false;
  }

  // INTERCHROMOSOMAL
  for(uint i=0; i!=movSize; i++) {
    double moveX=mov.mat(i,0);
    double moveY=mov.mat(i,1);
    double moveZ=mov.mat(i,2);
    double moveRadius=mov.chr[mov.start+i].getRadius();
    for(chriter = chromosomes.begin(); chriter != chromosomes.end(); chriter++) {
      if(*chriter != mov.chr) { 
        for(beaditer = chriter->begin(); beaditer != chriter->end(); beaditer++) {
          if(util::coordClash(beaditer->getX(),beaditer->getY(), beaditer->getZ(), moveX, moveY, moveZ, beaditer->getRadius()-0.00000001, moveRadius)) {
            return true;
          }
        }
      }
    }
  }
  return false;
}

bool Model::detectClash(){

  vector<Chromosome>::iterator chriter1; 
  vector<Chromosome>::iterator chriter2;  
  boost::container::static_vector<Bead,MAXSIZE>::iterator beaditer1;
  boost::container::static_vector<Bead,MAXSIZE>::iterator beaditer2;
  
  for(chriter1 = chromosomes.begin(); chriter1 != chromosomes.end(); chriter1++) {
    for(beaditer1 = chriter1->begin(); beaditer1 != chriter1->end(); beaditer1++) {
      for(chriter2 = chromosomes.begin(); chriter2 != chromosomes.end(); chriter2++) {
        for(beaditer2 = chriter2->begin(); beaditer2 != chriter2->end(); beaditer2++) {
          if(*beaditer1 != *beaditer2) { 
            
            if(util::coordClash(beaditer1->getX(),beaditer1->getY(), beaditer1->getZ(), beaditer2->getX(),beaditer2->getY(), beaditer2->getZ(), beaditer1->getRadius(), beaditer2->getRadius())) {              
	      return true;
            }
          }
        }
      }
    }
  }

  
  return false; 
}

uint weightedSamplePos(std::vector<double> weights, ENG &eng) {  
  std::vector<double> cumSumWeights(0, weights.size());
  double weightSum = accumulate(weights.begin(), weights.end(), 0);
  partial_sum(weights.begin(), weights.end(), weights.begin());
  double rnd = util::unif_real(0.0,weightSum, eng);
  uint pos = std::lower_bound(weights.begin(), weights.end(), rnd) - weights.begin();
  return pos;
}



Move Model::accept(Move &mov) { // Accepts the coordinates of mov, returns the coordinates of the altered beads (as a new Move).
  Move oldStructure(mov.chr);
  oldStructure.start = mov.start;
  oldStructure.stop = mov.stop;
  
  uint movSize = (mov.stop-mov.start)+1;
  CoordinateMatrix oldMat(movSize,3);
  for(uint i=0; i!=movSize; i++) {
    // Store old coordinates:
    oldMat(i,0) = mov.chr[mov.start+i].getX();
    oldMat(i,1) = mov.chr[mov.start+i].getY();
    oldMat(i,2) = mov.chr[mov.start+i].getZ();
    // Set new coordinates:
    mov.chr[mov.start+i].setCoordinates(mov.mat(i,0), mov.mat(i,1), mov.mat(i,2));
  }
  oldStructure.mat = oldMat;
  
  return oldStructure;
}


string Model::getCMM() {
  stringstream res;

  vector<Chromosome>::iterator chriter;
  boost::container::static_vector<Bead,MAXSIZE>::iterator beaditer;
  double smallestBead = numeric_limits<double>::infinity();
  for(chriter = chromosomes.begin(); chriter != chromosomes.end(); chriter++) {
    for(beaditer = chriter->begin(); beaditer != chriter->end(); beaditer++) {
      if(beaditer->getRadius() < smallestBead) {
	smallestBead = beaditer->getRadius();
      }
    }
  }
  
  double linkRadius = smallestBead * 0.5;

  res << "<marker_set name=\"" <<  this->name << "\">" << endl;
  uint counter = 0;
  string prevChr;
  uint prevId=0;
  for(chriter = chromosomes.begin(); chriter != chromosomes.end(); chriter++) {
    for(beaditer = chriter->begin(); beaditer != chriter->end(); beaditer++) {
      res << "<marker id=\"" << counter << "\" x=\"" << beaditer->getX() << "\" y=\"" << beaditer->getY() << "\" z=\"" << beaditer->getZ() << "\" radius=\"" << beaditer->getRadius() << "\" r=\"" << beaditer->color[0] << "\" g=\"" << beaditer->color[1] << "\" b=\"" << beaditer->color[2] << "\" chrID=\"" << chriter->getName() << "\" beadID=\"" << beaditer->getID() << "\"/>" << endl;

      if(prevChr == chriter->getName()) {
	res << "<link id1=\"" << prevId << "\" id2=\"" << counter << "\" r=\"" << beaditer->color[0] << "\" g=\"" << beaditer->color[1] << "\" b=\"" << beaditer->color[2] << "\" radius=\"" << linkRadius << "\"/>" << endl;
      }
      prevChr = chriter->getName();
      prevId = counter;
      counter++;
    }    
  }
  res << "</marker_set>" << endl;
  return res.str();
}

void Model::writeCMM(string fileName) {
  std::ofstream myfile;
  myfile.open(fileName.c_str());
  myfile << this->getCMM();
  myfile.close();
}

void Model::readGtrack(string filename, bool scaleBeadSizes/*=false*/, double nuclearOccupancy /*=0.2*/) {
  string line;
  ifstream myfile(filename.c_str());
  assert(myfile.is_open());

  vector<string> header;

  map<string,Chromosome> chromosomes;
  map<string, string> beadId2chr;
  vector<string> laminIds;
  vector<pair<string,string> > interactionIds;
  vector<pair<string,string> > interactionDistanceIds;
  vector<pair<string,string> > nonInteractionDistanceIds;
  std::map<string,double> distInfo;
  std::map<string,double> weightInfo;
  std::map<string,double> boundaryInfo;

  bool hasLamin = false;
  bool hasEdges = false;
  bool hasColor = false;
  
  while(getline(myfile,line)) {    
    if(line.substr(0,3) == "###") {
      header = util::split(line.substr(3),'\t');
      hasLamin = find(header.begin(), header.end(), "lamin") != header.end();
      hasEdges = find(header.begin(), header.end(), "edges") != header.end();
      hasColor = find(header.begin(), header.end(), "color") != header.end();      
    }
    if(line.substr(0,1) == "#") { // Ignore comments
      continue;
    }

    if(hasLamin) {
      assert(this->hasNucleus);
    }

    vector<string> dataLine = util::split(line,'\t');
    map<string, string> fieldParser = util::makeMap(header, dataLine);

    if (chromosomes.find(fieldParser["seqid"]) == chromosomes.end() ) {
      Chromosome chr(fieldParser["seqid"]);
      chromosomes[fieldParser["seqid"]] = chr;
    }

    Bead myBead(boost::lexical_cast<int>(fieldParser["start"]), boost::lexical_cast<int>(fieldParser["end"]), boost::lexical_cast<double>(fieldParser["radius"]), fieldParser["id"]);
    myBead.setCoordinates(0,0,0); // Coordinates are not defined at this point, set all to zero.
    
    beadId2chr[fieldParser["id"]] = fieldParser["seqid"];
    
    if(hasColor) {
      vector<string> col = util::split(fieldParser["color"], ',');
      myBead.setColor(boost::lexical_cast<float>(col[0]), boost::lexical_cast<float>(col[1]), boost::lexical_cast<float>(col[2]));
    }
    
    
    
    chromosomes[fieldParser["seqid"]].addBead(myBead);

    if(hasLamin and fieldParser["lamin"] == "1") {
      laminIds.push_back(fieldParser["id"]);
    }


    if(hasEdges and fieldParser["edges"] != ".") {
      vector<string> edgeList = util::split(fieldParser["edges"], ';');
      pair<string, string> idPair1;
      string b2id;
      for(uint i=0; i!=edgeList.size(); i++) {
        b2id = edgeList[i];
        if(b2id.find("=") == std::string::npos){ // if edge information does not contain edge weights 
          idPair1 = make_pair(fieldParser["id"], b2id);
          assert(idPair1.first != idPair1.second);	
 
          if(find(interactionIds.begin(), interactionIds.end(), idPair1) == interactionIds.end())  {
            interactionIds.push_back(idPair1); // idPair2 should not be added here, since addConstraint will add both pairs automatically
          }
        }
        else
        { // edge information contains edge weights
          b2id = edgeList[i].substr(0,edgeList[i].find("="));
          string edgeInfoAll =  edgeList[i].substr(edgeList[i].find("=")+1);
          idPair1 = make_pair(fieldParser["id"], b2id);
          assert(idPair1.first != idPair1.second);	
         
          if(edgeInfoAll.find(".") == 0){//No weighting on edges: treat like ordinary interaction
            
            if(find(interactionIds.begin(), interactionIds.end(), idPair1) == interactionIds.end())  {
              interactionIds.push_back(idPair1); // idPair2 should not be added here, since addConstraint will add both pairs automatically
            }
          }
          else{ // split edge information
            vector<double> edgeInfoDetail = util::splitDbl(edgeInfoAll,','); 
            assert(edgeInfoDetail.size()==3);
            assert(edgeInfoDetail[1] > 0);
            assert(edgeInfoDetail[2] == 0 or edgeInfoDetail[2]==1);
            // edgeInfoDetail contains [0]=distance (double), [1]=weight (<0,1]), [2]=boundary (0/1). 
            // If boundary is 1, define distance as the minimum required between the two "interacting" beads. (penalty for smaller distance than this defined distance)
            // If boundary is 0, use distance as the desired distance between the beads (penalty for larger distance than this defined distance)
            if(edgeInfoDetail[2]==0){
              interactionDistanceIds.push_back(idPair1);
            }
            else{
              nonInteractionDistanceIds.push_back(idPair1);
            } 
                     
            // Create bead_id pair first value as a string with comma as delimiter: "beadID1,beadID2"
            string myFirstValue = fieldParser["id"] + "," + b2id;
            distInfo[myFirstValue] = edgeInfoDetail[0];
            weightInfo[myFirstValue] = edgeInfoDetail[1];
            boundaryInfo[myFirstValue] = edgeInfoDetail[2];
          }
          
        }
    
      }
    }


  }
  myfile.close();
  

  // Adding the Chromosomes:
  for(map<string,Chromosome>::iterator it = chromosomes.begin(); it != chromosomes.end(); it++) {
    if(not hasColor) {
      double r = util::unif_real(0,1,randEngine);
      double g = util::unif_real(0,1,randEngine);
      double b = util::unif_real(0,1,randEngine);
      it->second.setColor(r,g,b);
    }
    this->addChromosome(it->second);
  }

  // Scale bead-sizes according to nuclear occupancy:
  if(scaleBeadSizes) {
    cerr << "# Rescaled bead sizes from " << this->getNuclearOccupancy();
    this->rescaleBeadSizes(nuclearOccupancy);
    cerr << " to " << this->getNuclearOccupancy() << endl;    
  }
  
  // Adding Interactions:
  for(vector<pair<string,string> >::iterator it = interactionIds.begin();  it != interactionIds.end(); it++) {
    //cout <<it->first << " " << it->second << endl;
    if(beadId2chr[it->first] == beadId2chr[it->second]) {
      this->addInteractionConstraint(it->first, it->second, 1, INTERACTION_INTRA);
    }
    else {
      this->addInteractionConstraint(it->first, it->second, 1, INTERACTION_INTER);
    }
  }

  // Adding Distance Interactions based on distances:
  for(vector<pair<string,string> >::iterator it = interactionDistanceIds.begin();  it != interactionDistanceIds.end(); it++) {
    
    // Obtaining distance and weight information
    string keyString = it->first + "," + it->second;
    std::map<string,double>::iterator itDist = distInfo.find(keyString);
    std::map<string,double>::iterator itWeight = weightInfo.find(keyString);
   
    
    this->addInteractionDistanceConstraint(it->first, it->second, itDist->second, itWeight->second, INTERACTION_DIST);
  }
  
    
  // Adding Forced Distance to NON interacting beads:
  for(vector<pair<string,string> >::iterator it = nonInteractionDistanceIds.begin();  it != nonInteractionDistanceIds.end(); it++) {
    
    // Obtaining distance and weight information
    string keyString = it->first + "," + it->second;
    std::map<string,double>::iterator itDist = distInfo.find(keyString);
    std::map<string,double>::iterator itWeight = weightInfo.find(keyString);
   
    
    this->addNonInteractionDistanceConstraint(it->first, it->second, itDist->second, itWeight->second, NON_INTERACTION_DIST);
  }

  // Adding lamin:
  
  for(vector<string>::iterator it = laminIds.begin();  it != laminIds.end(); it++) {
    this->addLaminConstraint(*it);
  }  

  // Initialize structures randomly:
  this->resetAllChromosomes(1000); //, true);

  uint N=this->getNumberOfBeads();
  cerr << "# beads: " << N << endl;
  cerr << "# interactions: " << interactionIds.size() << endl;
  cerr << "# interactions with given distance: " << interactionDistanceIds.size() << endl;
  cerr << "# non-interactions with given minimum distance: " << nonInteractionDistanceIds.size() << endl;
  cerr << "# lamin beads: " << laminIds.size() << endl;
  cerr << "# non-lamin beads: " << N-laminIds.size() << endl;  

}

double Model::getTotalChrLength() {
  vector<Chromosome>::iterator chriter;
  boost::container::static_vector<Bead,MAXSIZE>::iterator beaditer;
  double res = 0;
  for(chriter = chromosomes.begin(); chriter != chromosomes.end(); chriter++) {
    for(beaditer=chriter->begin(); beaditer!=chriter->end(); beaditer++) {
      res += beaditer->getRadius() * 2;
    }
  }
  return res;
}

uint Model::getNumberOfBeads() {
  vector<Chromosome>::iterator chriter;
  uint res = 0;
  for(chriter = chromosomes.begin(); chriter != chromosomes.end(); chriter++) {
    res += chriter->size();
  }
  return res;
}


uint Model::numberOfChromosomes() {
  return chromosomes.size();
}


void Model::resetAllChromosomes(uint maxAttempts /*=1000*/) { //, bool updateClash /*=false*/) {
  vector<Chromosome>::iterator chriter;
  boost::container::static_vector<Bead,MAXSIZE>::iterator beaditer;  
  double L = (this->getTotalChrLength() / this->numberOfChromosomes()) * 0.1;
  for(chriter = chromosomes.begin(); chriter != chromosomes.end(); chriter++) {
    vector<double> radii;
    Move mov(*chriter);
    boost::container::static_vector<Bead,MAXSIZE>::iterator beaditer;
    for(beaditer=chriter->begin(); beaditer!=chriter->end(); beaditer++) {
      radii.push_back(beaditer->getRadius());
    }
    mov.start = 0;
    mov.stop = chriter->size()-1;
    for(uint i=0; i!=maxAttempts; i++) {
      CoordinateVector origin(3);
      origin(0) = util::unif_real(-L, L, randEngine);
      origin(1) = util::unif_real(-L, L, randEngine);
      origin(2) = util::unif_real(-L, L, randEngine);
      CoordinateMatrix mat = util::createSelfAvoidingWalk(radii, origin, randEngine);
      
      if(mat.size1()!=0){
        mov.mat = mat;
        mov.validate();
        if (not isClashing(mov)) {
	  cerr << "# " << chriter->getName() << " done!" << endl;
	  Move old = this->accept(mov); //, updateClash);
	  break;
        }
      }
    }
  }  
}

void Model::addNucleusConstraints(double springConstant /*=1*/) {
  assert(this->hasNucleus);
  vector<Chromosome>::iterator chriter;
  boost::container::static_vector<Bead,MAXSIZE>::iterator beaditer;  
  for(chriter = chromosomes.begin(); chriter != chromosomes.end(); chriter++) {
    for(beaditer=chriter->begin(); beaditer!=chriter->end(); beaditer++) {
      this->addNucleusConstraint(beaditer->getID(),springConstant);
    }
  }
}


void Model::addBoundaryConstraints(double springConstant /*=1*/) {
  assert(this->hasNucleus);
  vector<Chromosome>::iterator chriter;
  boost::container::static_vector<Bead,MAXSIZE>::iterator beaditer;  
  for(chriter = chromosomes.begin(); chriter != chromosomes.end(); chriter++) {
    for(beaditer=chriter->begin(); beaditer!=chriter->end(); beaditer++) {
      this->addBoundaryConstraint(beaditer->getID(), this->nucleusRadius, springConstant);
    }
  }
}

void Model::addCenterConstraints(double springConstant /*=1*/) {
  assert(this->hasNucleus);
  vector<Chromosome>::iterator chriter;
  boost::container::static_vector<Bead,MAXSIZE>::iterator beaditer;  
  for(chriter = chromosomes.begin(); chriter != chromosomes.end(); chriter++) {
    for(beaditer=chriter->begin(); beaditer!=chriter->end(); beaditer++) {
      this->addCenterConstraint(beaditer->getID(), springConstant);
    }
  }
}



void Model::addSmartConstraints(double springConstant /*=1*/) {
  // Adds CenterConstraint to all non-lamin beads
  vector<Chromosome>::iterator chriter;
  boost::container::static_vector<Bead,MAXSIZE>::iterator beaditer; 
  map<string,vector<Constraint> >::iterator cmap_it;
  vector<Constraint>::iterator c_it;  
  bool haslamin;
  for(chriter = chromosomes.begin(); chriter != chromosomes.end(); chriter++) {
    for(beaditer=chriter->begin(); beaditer!=chriter->end(); beaditer++) {
      haslamin = false;
      cmap_it = constraints.find(beaditer->getID());      
      if(cmap_it != constraints.end()) { // If bead has constraints, loop through and look for 'LAMIN':
	for(c_it = cmap_it->second.begin(); c_it != cmap_it->second.end();c_it++) {
	  if(c_it->getConstraintType() == "LAMIN") {
	    haslamin = true;
	  }
	}
      }    
      if(!haslamin) {
        this->addCenterConstraint(beaditer->getID(),springConstant);
      }
    } 
  } 
}

std::vector<double> Model::calculateScorePerConstraint(Chromosome& chr, uint startPos, uint endPos){
  
  std::vector<double> myScore(C_T_COUNT,0.0);
  
  boost::container::static_vector<Bead,MAXSIZE>::iterator it;
  vector<Constraint>::iterator it2;
  for(it = chr.begin() + startPos; it != chr.begin() + endPos+1; it++) {
    for(it2 = constraints[it->getID()].begin(); it2 != constraints[it->getID()].end(); it2++) {
      if ((find(chr.begin()+startPos, chr.begin()+endPos+1, *(it2->bead1)) != chr.begin()+endPos+1) and (find(chr.begin()+startPos, chr.begin()+endPos+1, *(it2->bead2)) != chr.begin()+endPos+1) ) {
        myScore[it2->getConstraintTypeID()] += it2->eval()/2;
      }
      else {
        myScore[it2->getConstraintTypeID()] += it2->eval();
      }
    }
  }
  
  return myScore;
}

double Model::getTotalBeadVolume() {
  vector<Chromosome>::iterator chriter;
  boost::container::static_vector<Bead,MAXSIZE>::iterator beaditer;
  double res = 0;
  for(chriter = chromosomes.begin(); chriter != chromosomes.end(); chriter++) {
    for(beaditer=chriter->begin(); beaditer!=chriter->end(); beaditer++) {
      res += (4.0/3.0) * util::PI * pow(beaditer->getRadius(),3);
    }
  }
  return res;
}

double Model::getTotalGenomeLength() {
  vector<Chromosome>::iterator chriter;
  boost::container::static_vector<Bead,MAXSIZE>::iterator beaditer;
  double res = 0;
  for(chriter = chromosomes.begin(); chriter != chromosomes.end(); chriter++) {
    for(beaditer=chriter->begin(); beaditer!=chriter->end(); beaditer++) {
      res += beaditer->getBpSpan();
    }
  }
  return res;
}

void Model::rescaleBeadSizes(double occupancyFactor) {
  assert(this->hasNucleus);
  vector<Chromosome>::iterator chriter;
  boost::container::static_vector<Bead,MAXSIZE>::iterator beaditer;
  double beadVol = this->getTotalBeadVolume();
  double genomeLength = this->getTotalGenomeLength();
  for(chriter = chromosomes.begin(); chriter != chromosomes.end(); chriter++) {
    for(beaditer=chriter->begin(); beaditer!=chriter->end(); beaditer++) {
      double newRadius = pow((occupancyFactor * (beaditer->getBpSpan()) / (genomeLength)),1/3.0) * (this->nucleusRadius);
      beaditer->setRadius(newRadius);
    }
  }
}

double Model::getNuclearOccupancy() {
  assert(this->hasNucleus);
  double nuclearVolume = (4.0/3.0) * util::PI * pow(this->nucleusRadius,3);
  return this->getTotalBeadVolume() / nuclearVolume;
}
