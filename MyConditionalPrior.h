#ifndef DNest4_TestCase_MyConditionalPrior
#define DNest4_TestCase_MyConditionalPrior

#include "DNest4/code/DNest4.h"

// hyperparameters setting interim prior for object properties
class MyConditionalPrior:public DNest4::ConditionalPrior
{
  private:
    double f_lim; // flux limit
    double f_lim_min, f_lim_max;

    double gamma; // inverse slope for flux prior
    double gamma_min, gamma_max;

    // position limits
    double x_min, x_max;
    double y_min, y_max;

    double perturb_hyperparameters(DNest4::RNG& rng);

  public:
    MyConditionalPrior(
      double f_lim_min, double f_lim_max,
      double gamma_min, double gamma_max);

    void from_prior(DNest4::RNG& rng);

    double log_pdf(const std::vector<double>& vec) const;

    // implements the inverse CDF
    void from_uniform(std::vector<double>& vec) const;

    // implements the CDF
    void to_uniform(std::vector<double>& vec) const;

    void print(std::ostream& out) const;
};

#endif
