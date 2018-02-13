#ifndef DNest4_TestCase_Zernike
#define DNest4_TestCase_Zernike

#include <vector>
#include "DNest4/code/DNest4.h"

class Zernike
{
  private:
    int j_max; // maximum order

    std::vector<double> param; // parameter vector
    double param_min, param_max;

    // definition of radial polynomials
    std::vector<double> R_norm;
    std::vector< std::vector<double> > R_order;
    std::vector< std::vector<double> > R_coeff;

    // parameter indexing (OSA/ANSI)
    std::pair<int, int> index(int j) const;

  public:
    Zernike(int j_max, double param_min, double param_max);

    double evaluate(double x, double y) const;    

    void from_prior(DNest4::RNG& rng);
    double perturb(DNest4::RNG& rng);

    void print(std::ostream& out) const;
};

#endif
