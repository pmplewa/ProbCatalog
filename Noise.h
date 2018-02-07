#ifndef DNest4_TestCase_Noise
#define DNest4_TestCase_Noise

#include "../DNest4/code/DNest4.h"
#include "Zernike.h"

class UniformNoise
{
  private:
    double mean; // noise mean
    double mean_min, mean_max;

    double sigma;  // noise standard deviation
    double sigma_min, sigma_max;

  public:
    UniformNoise(
      double mean_min, double mean_max,
      double sigma_min, double sigma_max);

    void from_prior(DNest4::RNG& rng);
    double perturb(DNest4::RNG& rng);      

    double get_mean(double , double ) const { return mean; };
    double get_sigma() const { return sigma; };

    void print(std::ostream& out) const;
};

class VariableNoise
{
  private:    
    Zernike mean_function; // noise mean

    double sigma; // noise standard deviation
    double sigma_min, sigma_max;

  public:
    VariableNoise(
      double j_max, double param_min, double param_max,
      double sigma_min, double sigma_max);

    void from_prior(DNest4::RNG& rng);
    double perturb(DNest4::RNG& rng);    

    double get_mean(double x, double y) const { return mean_function.evaluate(x, y); };
    double get_sigma() const { return sigma; };

    void print(std::ostream& out) const;
};

#endif
