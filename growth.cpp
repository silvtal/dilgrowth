#include <RcppArmadilloExtensions/sample.h>
#include <Rcpp.h>
using namespace Rcpp ;
using namespace std ;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
NumericVector pick_new_bugs( NumericVector x, 
                              double size,
                              bool replace, 
                              NumericVector prob) {
  NumericVector pos = RcppArmadillo::sample(x, size, replace, prob);
  return (pos);
}


// [[Rcpp::export]]
NumericVector growth(NumericVector this_timestep, int abun_total, double grow_step) {
  
  // We will sample the growth positions from here
  NumericVector arr(this_timestep.size());
  for(int x = 0; x < this_timestep.size(); ++x)
    arr[x] = x;
  
  // Choose "step" (this is why we need "abun_total")
  double step;
  if ((sum(this_timestep) + grow_step) > abun_total) {
    // Avoid growing too much (when step>1)
    step = (abun_total - sum(this_timestep));
    
  } else if (sum(this_timestep) < grow_step) {
    // Ensure it is not too big of a step
    int half = trunc(sum(this_timestep)/2);
    step = std::max(half, 1);
    
  } else {
    // If grow_step is OK
    step = grow_step;
  }
  
  // Grow (loop: as many times as "step" indicates)
  NumericVector new_bugs = pick_new_bugs(arr, step, FALSE, this_timestep);
  int bug;
  for (std::size_t i = 0; i < new_bugs.size(); i++) {
    bug = new_bugs[i];
    this_timestep[bug] = (this_timestep[(bug)] + 1.0);
  }
  return(this_timestep);
}


// [[Rcpp::export]]
NumericVector full_growth(NumericVector this_timestep, int abun_total, double grow_step) {
  while (sum(this_timestep) < abun_total) {
    this_timestep = growth(this_timestep, abun_total, grow_step);
  }
  return(this_timestep);
}


// [[Rcpp::export]]
NumericVector growth_with_interactions(NumericVector this_timestep, int abun_total, double grow_step, NumericMatrix interactions) {
  
  // We will sample the growth positions from here
  NumericVector arr(this_timestep.size());
  for(int x = 0; x < this_timestep.size(); ++x)
    arr[x] = x;
  
  // Choose "step" (this is why we need "abun_total")
  double step;
  if ((sum(this_timestep) + grow_step) > abun_total) {
    // Avoid growing too much (when step>1)
    step = (abun_total - sum(this_timestep));
    
  } else if (sum(this_timestep) < grow_step) {
    // Ensure it is not too big of a step
    int half = trunc(sum(this_timestep)/2);
    step = std::max(half, 1);
    
  } else {
    // If grow_step is OK
    step = grow_step;
  }

  // Get final probability vector
  NumericVector prob = this_timestep/sum((this_timestep)); // abs abundances
  prob = wrap(as<arma::vec>(prob) + (as<arma::mat>(interactions) * as<arma::vec>(prob))); // x + A * x
  prob[prob < 0] = 0; // no negative probabilities
  
  // Grow (loop: as many times as "step" indicates)
  NumericVector new_bugs = pick_new_bugs(arr, step, FALSE, prob);
  int bug;
  for (std::size_t i = 0; i < new_bugs.size(); i++) {
    bug = new_bugs[i];
    this_timestep[bug] = (this_timestep[(bug)] + 1.0);
  }
  return(this_timestep);
}
