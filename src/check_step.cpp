#include <Rcpp.h>
using namespace Rcpp;

//' @export
// [[Rcpp::export]]
int check_step(NumericVector this_timestep,
                  int abun_total,
                  int grow_step) {
  // Choose "step" (this is why we need "abun_total")
  int step;
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
  return(step);
}
