#include <iostream>
#include "../src/Util.h"
#include "../src/Chromosome.h"
#include "../src/Bead.h"
#include "../src/Model.h"
#include "../src/Constraint.h"
#include "../src/Randomizer.h"
#include "../src/MCMC.h"
#include <assert.h>
#include <vector>
#include <math.h>
#include <limits>

#include <sstream>

#include "../src/tclap-1.2.1/include/tclap/CmdLine.h"

#define SSTR( x ) dynamic_cast< std::ostringstream & >( \
        ( std::ostringstream() << std::dec << x ) ).str()

using namespace std;

struct Args {
  string fileName;
  string outFile;
  uint seed;
  double radius;
  uint nIter;
  double maxtemp;
  double coolrate;
  double occupancy;
  uint verbose;
  string modelName;
  vector<string> excludeMoves;
  bool centerConstraint;
  bool nucleusConstraint;
  bool printStructures;
};

Args parseArguments(int argc, char** argv) {
  string VERSION="1.0.1";
  Args args;
    try {  
    TCLAP::CmdLine cmd("Chrom3D is a 3D genome modeling platform designed to incorporate a versatile set of constraints. For example, Chrom3D can simultaneously incorporate chromosomal interaction constraints (such as Hi-C) and constraints from chromosome association with the nuclear periphery." , ' ', VERSION);

    TCLAP::UnlabeledValueArg<string> nolabel( "filename", "Input file In gtrack format", true, "/dev/null", "Gtrack file"  );
    cmd.add( nolabel );	

    ////////
    /*
    ifstream f(nolabel.getValue().c_str());
    assert(f.good());
    */
    
    TCLAP::ValueArg<string> outfileArg("o","outstruct","Filename of the final, optimized structure in Chimera (CMM) format. If not specified, stdout will be used",false,"/dev/null","filename");

    TCLAP::ValueArg<uint> seedArg("s","seed","Seed used for randomization (default: 1234)",false,1234,"unsigned integer");

    TCLAP::ValueArg<double> radiusArg("r","radius","Radius of the nucleus in micrometer (default: 5.0)",false,5.0,"floating point number");

    TCLAP::ValueArg<uint> niterArg("n","iterations","Number of iterations (default: 2e+06)",false,2e+06,"unsigned integer");

    TCLAP::ValueArg<double> maxtempArg("T","maxtemp","Starting temperature of simulation (default: 1.0)",false,1.0,"floating point number");

    TCLAP::ValueArg<double> coolrateArg("c","coolrate","Cooling rate of temperature per iteration (default: 0)",false,0,"floating point number");

    TCLAP::ValueArg<double> occupancyArg("y","occupancy","Scale total volume of the model beads relative to the volume of the nucleus. Specified as a value between 0 and 1. (default: no scaling). Recommended value for e.g. human: 0.15.",false,0,"floating point number");

    TCLAP::ValueArg<uint> verboseArg("l","log","Write logging information to stderr every <number> iteration. This number must be between 0 and the total number of iterations (default: 0, not verbose)", false,0,"unsigned integer");

    TCLAP::ValueArg<string> modelnameArg("m","modelname","Name of model for the given run. (default: chrom3d_model)", false, "chrom3d_model" ,"string");

    vector<string> allowed;
    allowed.push_back("crankshaft");
    allowed.push_back("microcrank");
    allowed.push_back("armwiggle");
    allowed.push_back("armrotate");
    allowed.push_back("translate");
    allowed.push_back("rotate");

    TCLAP::ValuesConstraint<string> allowedVals( allowed );    
    TCLAP::MultiArg<string> excludeArg("e", "exclude", "Chromatin moves to exclude from the simulations", false, &allowedVals);


    //TCLAP::SwitchArg centerSwitch("","center","Add constraints such that all beads are pushed towards the center of the nucleus", false);
    TCLAP::SwitchArg nucleusSwitch("","nucleus","Add constraints such that all beads are pushed towards the inside of the nucleus", false);

    TCLAP::SwitchArg printfilesSwitch("","printmodels","Print intermediate models in CMM format (at intervals specified by -l/--log) with names consisting of the model name (specified by -m/--modelname) and the iteration number.", false);


    //////
    cmd.add( printfilesSwitch );
    cmd.add( excludeArg );
    //cmd.add( centerSwitch );
    cmd.add( nucleusSwitch );

    cmd.add( verboseArg );
    cmd.add( seedArg );
    cmd.add( occupancyArg );
    cmd.add( radiusArg );
    cmd.add( coolrateArg );
    cmd.add( maxtempArg );
    cmd.add( niterArg );
    cmd.add( modelnameArg );

    cmd.add( outfileArg );

    /////////////////
    
    cmd.parse( argc, argv );
    /// Add arguments to res:
    args.fileName = nolabel.getValue();
    args.outFile = outfileArg.getValue();
    args.seed = seedArg.getValue();
    args.radius = radiusArg.getValue();
    args.nIter = niterArg.getValue();
    args.maxtemp = maxtempArg.getValue();
    args.coolrate = coolrateArg.getValue();
    args.verbose = verboseArg.getValue();
    args.occupancy = occupancyArg.getValue();
    args.modelName = modelnameArg.getValue();

    //args.centerConstraint = centerSwitch.getValue();
    args.nucleusConstraint = nucleusSwitch.getValue();

      
    args.excludeMoves = excludeArg.getValue();

    args.printStructures = printfilesSwitch.getValue();
      
    assert(args.coolrate >= 0.0 and args.coolrate < 1.0);
    assert(args.occupancy >= 0.0 and args.occupancy < 1.0);
    
    assert(args.verbose >= 0 and args.verbose <= args.nIter);
    
  } catch (TCLAP::ArgException &e) { 
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; 
  }

    return(args);
}



int main(int argc, char** argv) {
  Args args = parseArguments(argc,argv);
  string fileName = args.fileName;
  
  Model model(args.modelName, args.seed);
  model.setNucleusRadius(args.radius);

  if(args.occupancy > 0) {
    model.readGtrack(fileName, true, args.occupancy);
  }
  else {
    model.readGtrack(fileName, false);
  }

  //if(args.centerConstraint) {
  //  model.addCenterConstraints();
  //}
  
  if(args.nucleusConstraint) {
    model.addNucleusConstraints();
  }

  
  if( find(args.excludeMoves.begin(), args.excludeMoves.end(), "crankshaft") ==  args.excludeMoves.end()) {
    model.addRandomizer(new CrankShaftRandomizer());
  }
  if( find(args.excludeMoves.begin(), args.excludeMoves.end(), "microcrank") ==  args.excludeMoves.end()) {
    model.addRandomizer(new MicroCrankShaftRandomizer());
  }
  if( find(args.excludeMoves.begin(), args.excludeMoves.end(), "armwiggle") ==  args.excludeMoves.end()) {
    model.addRandomizer(new ArmWiggleRandomizer());
  }
  if( find(args.excludeMoves.begin(), args.excludeMoves.end(), "armrotate") ==  args.excludeMoves.end()) {
    model.addRandomizer(new ArmRotateRandomizer());
  }
  if( find(args.excludeMoves.begin(), args.excludeMoves.end(), "translate") ==  args.excludeMoves.end()) {
    model.addRandomizer(new TranslateRandomizer(1.0));
  }
  if( find(args.excludeMoves.begin(), args.excludeMoves.end(), "rotate") ==  args.excludeMoves.end()) {
    model.addRandomizer(new RotationRandomizer());
  }
  
  //model.buildHashMap();

  MCMC mcmc(model);
  bool success;
  uint Niter=args.nIter;
  double maxTemp = args.maxtemp;
  double coolRate = args.coolrate; // Number between 0 and 1, giving the rate of cooling at each step
  double temp=maxTemp;

  uint verbosity = (args.verbose == 0) ? args.nIter+1 : args.verbose;

   
  if(args.verbose != 0) {
    cerr << "0 " << model.getLossScore(INTERACTION_INTRA) << " " << model.getLossScore(INTERACTION_INTER) << " " << model.getLossScore(PERIPHERY) << " " << model.getLossScore(CENTER) << " " << model.getLossScore(BOUNDARY) << " " << model.getLossScore(NUCLEUS) << " " << model.getLossScore(INTERACTION_DIST) << " " << model.getLossScore() << endl;

  }
  
  
  for(uint i=1; i!=args.nIter+1; i++) {
    temp *= (1-coolRate);
    success = mcmc.doMetropolisHastingsStep(temp); //simulated annealing when T < 1. When T = 1 that is Metropolis Hastings
    assert(success);
    if (i % verbosity == 0) {
      if(args.printStructures) {
	model.writeCMM(args.modelName + "_iter_" + SSTR(i) + ".cmm");
      }
      cerr << i << " " << model.getLossScore(INTERACTION_INTRA) << " " << model.getLossScore(INTERACTION_INTER) << " " << model.getLossScore(PERIPHERY) << " " << model.getLossScore(CENTER) << " " << model.getLossScore(INTERACTION_DIST) << " " << model.getLossScore() << endl;
    }
  }

  if(args.outFile != "/dev/null") {
    model.writeCMM(args.outFile);
  }
  else {
    cout << model.getCMM();
  }

  exit(EXIT_SUCCESS);
}
