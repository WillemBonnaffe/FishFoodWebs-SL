//* MCMCMH ALGORITHM *//

// goal: code an MCMC algorithm implemented with MH

// updates:
// @ v0.1 - 17/11/2019 - implemented automatic cacluation of gamma
// @ v0.2 - 25/01/2024 - implemented using pre-existing chain as starting point 


/* ********** */
/* INITIATION */
/* ********** */

// C++ dependencies
#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
List DEMCpp(List argList)
{
  // Unwrap arguments
  Function dTarget = as<Function>(argList["dTarget"]);
  float epsilon = as<float>(argList["epsilon"]);
  int nIt = as<int>(argList["nIt"]);

  // Initialise primer chain
  NumericMatrix primerChain = as<NumericMatrix>(argList["primerChain"]);
  int nPrimer = primerChain.nrow();
  int d = primerChain.ncol() - 1;

  // Initialize chain with primer
  NumericMatrix chain(nIt + nPrimer, d + 1);
  for (int k = 0; k < nPrimer; k++)
  {
    NumericVector Theta_local = primerChain(k, _); Theta_local.erase(0);
    float target_local = primerChain(k, 0);
    NumericVector chain_local = Theta_local; chain_local.push_front(target_local);
    chain(k,_) = chain_local;
  }

  // Initialise local variables 
  NumericVector Theta_local = chain(nPrimer - 1, _); Theta_local.erase(0);
  float target_local = chain(nPrimer - 1, 0);
  NumericVector chain_local = Theta_local; chain_local.push_front(target_local);

  // Set step size
  float gamma = 2.38 / pow(2.0 * Theta_local.size(), 0.5);

  // Iterations
  float r;
  int accepted = 0;
  int count = 0;

  for (int k = nPrimer; k < nIt + nPrimer; k++)
  {
    // Set distance vector with lag
    NumericVector delta = chain(as<int>(runif(1, k / 2, k)), _) - chain(as<int>(runif(1, k / 2, k)), _);
    delta.erase(0);

    // Update theta
    NumericVector Theta_newLocal = Theta_local + gamma * delta + epsilon * as<NumericVector>(runif(d, -1.0, 1.0));

    // Update target
    float target_newLocal = as<float>(dTarget(Theta_newLocal));

    // Metropolis ratio
    r = exp(target_newLocal - target_local);

    // Metropolis test
    if (drand48() < r)
    {
      accepted += 1;
      target_local = target_newLocal;
      Theta_local = Theta_newLocal;
    }

    // Update the chain
    chain_local = Theta_local; chain_local.push_front(target_local);
    chain(k,_) = chain_local;

    // Message
    if (count == (int)((nIt + nPrimer)/ 10))
    {
      Rcout << k << "/" << nIt + nPrimer << " | " << chain(k, 0) << " | " << chain(k, 1) << " | " << chain(k, 2) << "\n";
      count = 0;
    }
    count += 1;
  }

  // Return results
  Rcout << "p = " << (float)accepted / (float)(nIt + nPrimer) << "\n";
  return List::create(
    Named("dTarget") = dTarget(Theta_local),
    Named("Theta_0") = Theta_local,
    Named("gamma") = gamma,
    Named("chainList") = chain
  );
}
