#include <RcppArmadilloExtensions/sample.h>
#include <Rcpp/stats/random/rbinom.h>
#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// [[Rcpp::depends(RcppArmadillo)]]
//' @export
//' Sampling function for growth()
// [[Rcpp::export]]
NumericVector pick_new_bugs(NumericVector x,
                             double size,
                            bool replace,
                            NumericVector prob) {
  NumericVector pos = RcppArmadillo::sample(x, size, replace, prob);
  return (pos);
}

//' @export
//' This Rcpp function simulates the growth of a population of organisms over a
//' specified time step. It takes the current population abundance, the maximum
//' growth rate, and the maximum total abundance as inputs, and returns an
//' updated abundance vector that reflects the growth of the population. The
//' function calculates the growth rate based on the current population
//' abundance and the maximum possible total abundance, and calculates
//' the probability of growth for each organism in the population based on their
//' current abundance and any specified interactions between organisms. The
//' function then randomly selects organisms to grow based on these probabilities
//' and increases their abundance by grow_step (1 by default).
//' The function takes four arguments:
//' @param this_timestep A numeric vector representing the current abundance of
//' each organism in the population.
//' @param abun_total An integer representing the maximum total abundance that
//' the population can reach.
//' @param grow_step An integer representing the maximum growth rate of the
//' population, 1 by default.
//' @param interactions An optional numeric matrix representing the interaction
//' between organisms in the population. This argument is set to R_NilValue by default.
// [[Rcpp::export]]
NumericVector growth(NumericVector this_timestep,
                     int grow_step,
                     Rcpp::Nullable<Rcpp::NumericMatrix> interactions = R_NilValue) {

  // We will sample the growth positions from here
  NumericVector arr(this_timestep.size());
  for(int x = 0; x < this_timestep.size(); ++x)
    arr[x] = x;

  // Account for interactions if present
  NumericVector prob = this_timestep/sum((this_timestep)); // abs abundances
  if (interactions.isNotNull()) {
    prob = wrap(as<arma::vec>(prob) + (as<arma::mat>(interactions) * as<arma::vec>(prob))); // x * (r + A * x) == rx + A·x*2
    // 26x1                      // 26x26                  // 26x1  // diag of A here is 0. Summing here is equivalent to having a diagonal of 1.
  } else {
    prob = as<arma::vec>(prob); // x * (r + A * x) == rx + A·x*2
  }
  prob[prob < 0] = 0; // no negative probabilities

  // Grow (loop: as many times as "step" indicates)
  NumericVector new_bugs = pick_new_bugs(arr, grow_step, FALSE, prob);
  int bug;
  for (std::size_t i = 0; i < new_bugs.size(); i++) {
    bug = new_bugs[i];
    this_timestep[bug] = (this_timestep[(bug)] + 1.0);
  }
  return(this_timestep);
}


//' @export
//' This function consists in a loop that runs growth() as many times as needed
//' to reach a given population size.
// [[Rcpp::export]]
NumericVector full_growth(NumericVector this_timestep, int abun_total, int grow_step) {
  while (sum(this_timestep) < abun_total) {
    grow_step = check_step(this_timestep, abun_total, grow_step);
    this_timestep = growth(this_timestep, grow_step);
  }
  return(this_timestep);
}



// [[Rcpp::depends(RcppArmadillo)]]
//' This function simulates growth in a community by looking at the carrying
//' capacities of the group they belong to. It takes a named vector,
//' carrying_capacities.
//'
//' Growth is RANDOM: all species have the opportunity to grow each
//' time this function is called, but whether they will is random and depends on
//' their abundance (random). They grow at different rates depending not
//' only on which group they belong to but also on their carrying capacity.
//'
//' As opposed to growth(), the total growth rate is defined by the sum of the
//' growth rates of every group, which are defined as
//' \[grow_step * group's % of total carrying capacity]. Every group will have
//' grow_step (by default 1) of its members grow in each run, and they can be
//' from the same species or not. Also, how much will each species grow is
//' proportional to the carrying capacity of its group. This is to avoid group
//' extinction and also to ensure growth has a similar, proportional rate for
//' each group so all groups reach their CC at the same time.
//'
//' @export
// [[Rcpp::export]]
NumericVector growth_per_group(NumericVector x,
                               NumericVector carrying_capacities,
                               int grow_step,
                               Rcpp::Nullable<Rcpp::NumericMatrix> interactions = R_NilValue) { // named vector
  // Set variables
  double size = carrying_capacities.length();
  NumericVector newx = x;
  NumericVector prob = x/sum((x));
  CharacterVector names  = carrying_capacities.attr("names");
  CharacterVector groups = unique(names);

  // Account for interactions if present
  if (interactions.isNotNull()) {
    prob = wrap(as<arma::vec>(prob) + (as<arma::mat>(interactions) * as<arma::vec>(prob))); // x * (r + A * x) == rx + A·x*2
    // 26x1                      // 26x26                  // 26x1  // diag of A here is 0. Summing here is equivalent to having a diagonal of 1.
  } else {
    prob = as<arma::vec>(prob); // x * (r + A * x) == rx + A·x*2
  }

  // Growth per group
  for (String group : groups) {
    // indexing vector for this group
    LogicalVector g(size); // size is size of x! groups.size() is the # of groups
    for (int i = 0; i < size; i++) {
      g[i] = (names[i] == group);
    }
   // (only bugs from this group can grow)
    NumericVector group_prob(size);
    for (int i = 0; i < size; i++) {
      if (g[i]) {
        group_prob[i] = prob[i];
      } else {
        group_prob[i] = 0;
      }
    }
    NumericVector new_bugs = pick_new_bugs(x, grow_step, FALSE, group_prob);
    int bug;
    for (std::size_t i = 0; i < new_bugs.size(); i++) {  // grow_step times
      bug = new_bugs[i];
      x[bug] = (x[bug] + (1.0 * carrying_capacities[LogicalVector (carrying_capacities.attr("names")==groups[i])][0])); // % of total CC
      }
    }
  return newx;
}



// [[Rcpp::depends(RcppArmadillo)]]
//' This function simulates growth in a community by looking at the carrying
//' capacities of the group they belong to. It takes a named vector,
//' carrying_capacities.
//'
//' Growth is RANDOM but logistic: all species have the opportunity to grow each
//' time this function is called, but whether they will is random and depends on
//' their abundance (random). They grow at different rates depending not
//' only on which group they belong to but also to how close they are to their
//' carrying capacity (logistic)
//'
//' As opposed to growth(), the growth rate is given by a logistic function and
//' not grow_step. Another difference is that growing species are not chosen one
//' by one by sampling, but with a binomial function. That means that the number
//' of different species that can grow in each iteration of this function is not
//' limited.
//'
//' The difference with growth_by_group() is that growth is logistic in this
//' function.
//'
//' If the carrying capacity for a group was surpassed before starting the
//' growth cycle, the species of that group will die at a proportionate rate,
//' while the others may grow.
//'
//' @export
// [[Rcpp::export]]
NumericVector growth_log(NumericVector x,
                         NumericVector carrying_capacities,
                         Rcpp::Nullable<Rcpp::NumericMatrix> interactions = R_NilValue) { // named vector
  // Set variables
  double size = carrying_capacities.length();
  NumericVector newx = x;
  NumericVector prob = x/sum((x));
  CharacterVector names  = carrying_capacities.attr("names");
  CharacterVector groups = unique(names);

  // Get new abundance for each bug, group by group.
  // This is were we apply the probabilities: within each group, each species
  // might or might not grow.
  for (String group : groups) {
    // indexing vector for this group
    LogicalVector g(size);
    for (int i = 0; i < size; i++) {
      g[i] = (names[i] == group);
    }
    // sum_group => total abundance of the group
    double sum_group = 0.0;
    for (int i = 0; i < size; i++) {
      if (g[i]) {
        sum_group += x[i];
      }
    }

    // Account for interactions if present
    if (interactions.isNotNull()) {
      prob = wrap(as<arma::vec>(prob) + (as<arma::mat>(interactions) * as<arma::vec>(prob))); // x * (r + A * x) == rx + A·x*2
      // 26x1                      // 26x26                  // 26x1  // diag of A here is 0. Summing here is equivalent to having a diagonal of 1.
    } else {
      prob = as<arma::vec>(prob); // x * (r + A * x) == rx + A·x*2
    }

    // Here, we randomly determine whether each species in the group will grow
    for (int i = 0; i < size; i++) {
      if (g[i]) {
        if (rbinom(1, 1, prob[i])[0] == 1) { // more likely to grow or die if more abundant <- nothing if 0
          double growth_step(2 * (1 - sum_group / carrying_capacities[i]));
          newx[i] = x[i] + max(NumericVector(1, growth_step));
        }
      }
    }
  }
  return newx;
}
