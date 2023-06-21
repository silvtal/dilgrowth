#include <RcppArmadilloExtensions/sample.h>
#include <Rcpp/stats/random/rbinom.h>
#include <Rcpp.h>
using namespace Rcpp;
using namespace std;


// [[Rcpp::depends(RcppArmadillo)]]
//' Sampling function for growth_one_group(), growth()
//' @param arr is a vector of positions. If the abundance vector of interest has
//' 5 members, x must be [1 2 3 4 5].
//' @export
// [[Rcpp::export]]
NumericVector pick_new_bugs(NumericVector arr,
                            int size,
                            bool replace,
                            NumericVector prob) {
  NumericVector pos = RcppArmadillo::sample(arr, size, replace, prob);
  return (pos);
}

//' This Rcpp function simulates the growth of a community of organisms over a
//' specified time step. It takes the current community abundances, a fixed
//' growth step, and an optional interactions matrix as inputs, and returns
//' an updated abundances vector that reflects the growth of the community. The
//' function calculates the probability of growth for each organism based on
//' their current abundance and any specified interactions between organisms. The
//' function then randomly selects organisms to grow based on these probabilities
//' and increases their abundance by grow_step (1 by default).
//' The function takes four arguments:
//' @param this_timestep A numeric vector representing the current abundance of
//' each organism in the population.
//' @param grow_step An integer representing the growth step for organisms in
//' the community, 1 by default.
//' @param interactions An optional numeric matrix representing the interaction
//' between organisms in the population. This argument is set to R_NilValue by
//' default.
//' @export
// [[Rcpp::export]]
NumericVector growth_one_group(NumericVector this_timestep,
                               double grow_step,
                               Rcpp::Nullable<Rcpp::NumericMatrix> interactions = R_NilValue) {

  // We will sample the growth positions from here
  NumericVector arr(this_timestep.size());
  for(int x = 0; x < this_timestep.size(); ++x)
    arr[x] = x;

  // Account for interactions if present
  NumericVector prob = this_timestep/sum((this_timestep)); // abs abundances
  if (interactions.isNotNull()) {
    prob = wrap(as<arma::vec>(prob) + (as<arma::mat>(interactions.get()) * as<arma::vec>(prob)));
    for (int i = 0; i < prob.size(); i++) { // no negative probabilities allowed
      if (prob[i] < 0) {
        prob[i] = 0;
      }
    }
  }

  // Grow (loop: as many times as "step" indicates)
  NumericVector new_bugs = pick_new_bugs(arr, grow_step, true, prob);
  int bug;
  for (std::size_t i = 0; i < new_bugs.size(); i++) {
    bug = new_bugs[i];
    this_timestep[bug] = (this_timestep[(bug)] + 1.0);
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
//' their abundance. Their growth probability varies depending not only on their
//' abundances but also on their carrying capacity. Growth is not logistic.
//'
//' As opposed to growth_one_group(), the growth rate (probability of being
//' picked for growth) is determined for each individual by the carrying
//' capacity of its group(\[grow_step * group's % of total carrying capacity]).
//' An optionally passed interactions matrix has an effect as well. Every group
//' will have grow_step (by default 1) of its present individuals grow in each
//' run, and they can be from the same ASV/OTU/species or not. Also, growth step
//' is proportional to the carrying capacity of its group. This is to avoid group
//' extinction and also to ensure growth has a similar, proportional rate for
//' each group so all groups reach their CC at the same time.
//'`
//' @export
// [[Rcpp::export]]
NumericVector growth(NumericVector x,
                     NumericVector carrying_capacities,
                     int grow_step,
                     Rcpp::Nullable<Rcpp::NumericMatrix> interactions = R_NilValue) { // named vector

  // Set variables
  double size = carrying_capacities.length();
  NumericVector prob = x/sum((x));
  CharacterVector names  = carrying_capacities.attr("names");
  CharacterVector groups = unique(names);
  NumericVector group_ccs = carrying_capacities[groups];
  NumericVector arr(x.size());
  for (int i = 0; i < x.size(); ++i) {
    arr[i] = i;
  }

  // Account for interactions if present
  if (interactions.isNotNull()) {
    prob = wrap(as<arma::vec>(prob) + (as<arma::mat>(interactions.get()) * as<arma::vec>(prob)));
    for (int i = 0; i < prob.size(); i++) { // no negative probabilities allowed
      if (prob[i] < 0) {
        prob[i] = 0;
      }
    }
  }

  // Stop if all probabilities become 0
  if (sum(prob)==0) {
    Rf_error("After accounting for interactions, all probabilities became 0. No growth possible. Check the interaction matrix.");
  }

  // Growth per group
  for (String group : groups) {
    // indexing vector for this group
    LogicalVector g(size); // size is size of x! groups.size() is the # of groups
    for (int i = 0; i < size; i++) {
      g[i] = (names[i] == group);
    }
   // (only bugs from this group can grow)
   // (also check magnitude of grow_step)
    NumericVector group_prob(size);
    double sum_group = 0;
    for (int i = 0; i < size; i++) {
      if (g[i]) {
        group_prob[i] = prob[i];
        sum_group += x[i];
      } else {
        group_prob[i] = 0;
      }
    }
    sum_group = trunc(sum_group);
    int step = min(static_cast<int>(grow_step), static_cast<int>(sum_group)); // TODO

    // Stop if no one can grow in this group
    if (sum(group_prob)==0) {
      continue;
    }

    // Grow
    NumericVector new_bugs = pick_new_bugs(arr, step, true, group_prob);
    int bug;
    for (int i = 0; i < new_bugs.size(); i++) {  // grow_step times
      bug = new_bugs[i];
      x[bug] += (1.0 * (group_ccs[group]/sum(group_ccs))); // % of total CC
      }
    }
  return x;
}



//' This function simulates growth in a community by looking at the carrying
//' capacities of the group they belong to. It takes a named vector,
//' carrying_capacities.
//'
//' Growth is RANDOM but logistic: all species have the opportunity to grow each
//' time this function is called, but whether they will is random and depends on
//' their abundance (random). Their growth probability varies depending not
//' only on their abundances but also on their carrying capacity.
//'
//' As opposed to growth_one_group(), the growth rate is given by a logistic
//' function and not grow_step. Another difference is that growing species are
//' not chosen one by one by sampling, but with a binomial function. That means
//' that the number of different species that can grow in each iteration of this
//' function is not limited.
//'
//' The difference with growth() is that growth is logistic in this function and
//' there is not grow_step here.
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
      prob = wrap(as<arma::vec>(prob) + (as<arma::mat>(interactions.get()) * as<arma::vec>(prob)));
      for (int i = 0; i < prob.size(); i++) { // no negative probabilities allowed
        if (prob[i] < 0) {
          prob[i] = 0;
        }
      }
    }

    // Here, we randomly determine whether each species in the group will grow
    for (int i = 0; i < size; i++) {
      if (g[i]) {
        if (rbinom(1, 1, prob[i])[0] == 1) { // more likely to grow or die if more abundant <- nothing if 0
          double growth_step(2 * (1 - sum_group / carrying_capacities[i]));
          newx[i] += max(NumericVector(1, growth_step));
        }
      }
    }
  }
  return newx;
}


//' Checks if a given grow_step is ok for running a growth() function and adjusts
//' it accordingly if it's not (for example not allowing for it to cause too big
//' of a growth, surpassing total wanted final community abundance).
//' @param is_grow_step_a_perc Boolean: if false, grow_step is taken as a fixed
//' value, so the step will always be the same. If true, it is taken to indicate
//' a percentage - the step will be changed proportionally to the community size.
//' If grow_step is 0.02, 2% of the members in the community will grow the next
//' iteration
//' @export
// [[Rcpp::export]]
int check_step(NumericVector this_timestep,
               int abun_total,
               double grow_step,
               bool is_grow_step_a_perc = false) {

  int step = grow_step;
  if (is_grow_step_a_perc) {
    // Apply a percentage first if necessary
    if ((grow_step > 1) || (grow_step < 0)) {
      Rf_error("If grow_step is a percentage, it has to be a value between 0 (0% of the members of the community will grow) and 1 (100% of members will grow, the total abundance will duplicate every iteration).");
    } else {
      step = trunc(grow_step * sum(this_timestep));
      step = max(1, step); // never less than 1 !!
    }
  }

  // Check magnitude of step
  if (trunc((sum(this_timestep)) + grow_step) > abun_total) { // trunc: this_timestep might not have whole numbers.
    // Avoid growing too much (when step>1)
    step = (abun_total - sum(this_timestep));
  } else if (sum(this_timestep) < grow_step) {
    // Ensure there are enough bugs to grow with that step
    int half = trunc(sum(this_timestep)/2);
    step = std::max(half, 1);
  }
  return(step);
}


//' This function consists in a loop that runs growth(...)() functions as many
//' times as needed to reach a given population size.
//' @export
// [[Rcpp::export]]
NumericVector full_growth(NumericVector this_timestep,
                          int abun_total,
                          int grow_step,
                          bool is_grow_step_a_perc = false,
                          String func = "growth",
                          Rcpp::Nullable<Rcpp::NumericMatrix> interactions = R_NilValue,
                          Rcpp::Nullable<Rcpp::NumericVector> carrying_capacities = R_NilValue) {
  if (func == "growth_one_group") {
    while (sum(this_timestep) < abun_total) {
      int step      = check_step(this_timestep, abun_total, grow_step, is_grow_step_a_perc);
      this_timestep = growth_one_group(this_timestep, step, interactions.get());
    }
  } else if (func == "growth" && carrying_capacities.isNotNull()) {
    while (sum(this_timestep) < abun_total) {
      int step      = check_step(this_timestep, abun_total, grow_step, is_grow_step_a_perc);
      this_timestep = growth(this_timestep, carrying_capacities.get(), step, interactions.get());
    }
  } else if (func == "growth_log" && carrying_capacities.isNotNull()) {
    while (round(sum(this_timestep)) < abun_total) {
      this_timestep = growth_log(this_timestep, carrying_capacities.get(), interactions.get());
    }
  }
  return(this_timestep);
}
