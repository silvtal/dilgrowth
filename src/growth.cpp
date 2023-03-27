#include <RcppArmadilloExtensions/sample.h>
#include <Rcpp/stats/random/rbinom.h>
#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// [[Rcpp::depends(RcppArmadillo)]]
//' @export
// [[Rcpp::export]]
NumericVector pick_new_bugs( NumericVector x,
                              double size,
                              bool replace,
                              NumericVector prob) {
  NumericVector pos = RcppArmadillo::sample(x, size, replace, prob);
  return (pos);
}

//' @export
// [[Rcpp::export]]
NumericVector growth(NumericVector this_timestep,
                     int abun_total,
                     double grow_step,
                     Rcpp::Nullable<Rcpp::NumericMatrix> interactions = R_NilValue) {

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
  if (interactions.isNotNull()) {
    prob = wrap(as<arma::vec>(prob) + (as<arma::mat>(interactions) * as<arma::vec>(prob))); // x * (r + A * x) == rx + A·x*2
    // 26x1                      // 26x26                  // 26x1  // diag of A here is 0. Summing here is equivalent to having a diagonal of 1.
  } else {
    prob = as<arma::vec>(prob); // x * (r + A * x) == rx + A·x*2
  }
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


//' @export
// [[Rcpp::export]]
NumericVector full_growth(NumericVector this_timestep, int abun_total, double grow_step) {
  while (sum(this_timestep) < abun_total) {
    this_timestep = growth(this_timestep, abun_total, grow_step);
  }
  return(this_timestep);
}







// [[Rcpp::depends(RcppArmadillo)]]
//' This function simulates growth in a community by looking at the carrying
//' capacities of each species, which they have in common with other species in
//' their group. It takes a named vector, carrying_capacities.
//'
//' Growth is RANDOM but logistic: all species have the opportunity to grow each
//' time this function is called, but whether they will is random and depends on
//' their abundance (random), but they grow at different rates depending not
//' only on which group they belong to but also to how close they are to their
//' carrying capacity  (logistic)
//'
//' If the carrying capacity for a group was surpassed already, the species of
//' that group will die at a proportionate rate, while the others may grow.
//'
//' @export
// [[Rcpp::export]]
NumericVector growth_log(NumericVector x,
                         NumericVector carrying_capacities) { // named vector
  // Set variables
  double size = carrying_capacities.length();
  NumericVector newx(size);
  NumericVector prob = x/sum((x));
  CharacterVector names  = carrying_capacities.attr("names");
  CharacterVector groups = unique(names);

  // Get new abundance for each bug, group by group.
  // This is were we apply the probabilities: within each group, each species
  // might or might not grow.
  for (String group : groups) {
    // indexing vector for this group
    int n = names.size();
    LogicalVector g(size);
    for (int i = 0; i < n; i++) {
      g[i] = (names[i] == group);
    }
    // sum_group => total abundance of the group
    double sum_group = 0.0;
    for (int i = 0; i < size; i++) {
      if (g[i]) {
        sum_group += x[i];
      }
    }
    // Here, we randomly determine whether each species in the group will grow
    for (int i = 0; i < size; i++) {
      if (g[i]) {
        if (rbinom(1, 1, prob[i])[0] == 1) { // more likely to grow or die if more abundant <- nothing if 0
          NumericVector growth_step(1, (2 * (1 - sum_group / carrying_capacities[i])));
          newx[i] = x[i] + max(growth_step);
            ;
        } else {
          newx[i] = x[i];
        }
      }
    }
  }
  return newx;
}
