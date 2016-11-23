#ifndef GUARD_MCMC_H
#define GUARD_MCMC_H
#include "Model.h"
#include <math.h>
#include "Util.h"
#include <fstream>
#include <iostream>

class MCMC {
 public:
  MCMC(Model& mod);
  Model &model;
  bool doMetropolisHastingsStep(double T=1, uint maxAttempts=10000);
};

#endif
