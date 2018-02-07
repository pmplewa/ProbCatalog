#include "Noise.h"

using namespace std;

UniformNoise::UniformNoise(
  double mean_min, double mean_max,
  double sigma_min, double sigma_max)
:mean_min(mean_min)
,mean_max(mean_max)
,sigma_min(sigma_min)
,sigma_max(sigma_max)
{

}

void UniformNoise::from_prior(DNest4::RNG& rng)
{
  mean = mean_min + (mean_max-mean_min)*rng.rand();
  sigma = exp(log(sigma_min) + log(sigma_max/sigma_min)*rng.rand());
}

double UniformNoise::perturb(DNest4::RNG& rng)
{
  double logH = 0.;

  int which = rng.rand_int(2);

  if(which == 0)
  {
    // perturb noise mean
    mean += (mean_max-mean_min)*rng.randh();
    DNest4::wrap(mean, mean_min, mean_max);
    // logH += 0;  
  }
  else if(which == 1)
  {
    // perturb noise standard deviation
    double log_sigma = log(sigma);
    log_sigma += log(sigma_max/sigma_min)*rng.randh();
    DNest4::wrap(log_sigma, log(sigma_min), log(sigma_max));
    sigma = exp(log_sigma);
    // logH += 0;
  }

  return logH;
}

void UniformNoise::print(std::ostream& out) const
{
  out << mean << " " << sigma << " ";
}

VariableNoise::VariableNoise(
  double j_max, double param_min, double param_max,
  double sigma_min, double sigma_max)
:mean_function(Zernike(j_max, param_min, param_max))
,sigma_min(sigma_min)
,sigma_max(sigma_max)
{

}

void VariableNoise::from_prior(DNest4::RNG& rng)
{
  mean_function.from_prior(rng);
  sigma = exp(log(sigma_min) + log(sigma_max/sigma_min)*rng.rand());
}

double VariableNoise::perturb(DNest4::RNG& rng)
{
  double logH = 0.;

  int which = rng.rand_int(2);

  if(which == 0)
  {
    // perturb noise mean
    logH += mean_function.perturb(rng);
  }
  else if(which == 1)
  {
    // perturb noise standard deviation
    double log_sigma = log(sigma);
    log_sigma += log(sigma_max/sigma_min)*rng.randh();
    DNest4::wrap(log_sigma, log(sigma_min), log(sigma_max));
    sigma = exp(log_sigma);
    // logH += 0;
  }

  return logH;
}

void VariableNoise::print(std::ostream& out) const
{
  mean_function.print(out);
  out << sigma << " ";
}
