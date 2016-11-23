#include "MCMC.h"


MCMC::MCMC(Model& mod): model(mod) { }
	   
bool MCMC::doMetropolisHastingsStep(double T /*=1*/, uint maxAttempts /*=10000*/) {
  uint i=0;
  bool success = false;
  while(i < maxAttempts and !success) {
    try {
      Move proposedMove = model.proposeMove(); // Idea: if this fails (no move proposed), we could set the proposed Move to the same (identity) move.
      double before = model.getPartialLossScore(proposedMove.chr, proposedMove.start, proposedMove.stop); // score before move
      Move oldCoords = model.accept(proposedMove); //, true);
      double after = model.getPartialLossScore(proposedMove.chr, proposedMove.start, proposedMove.stop); // score after move
      if(after > before) {
	double pAccept = exp((before-after)/T);
	bool accept = util::unif_real(0,1, model.randEngine) < pAccept;
	if(not accept) {
	  Move tmp = model.accept(oldCoords); //, true); // roll back to previous coordinates.
	}
      }
      success = true;
    }
    catch(util::MoveException& e){
      i++;
    }
  }
  return success;
}

