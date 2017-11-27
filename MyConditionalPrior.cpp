#include "MyConditionalPrior.h"
#include "DNest4/code/DNest4.h"
#include "Data.h"
#include "Distributions.h"
#include <cmath>

using namespace std;

MyConditionalPrior::MyConditionalPrior(
  double f_lim_min, double f_lim_max,
  double gamma_min, double gamma_max)
:f_lim_min(f_lim_min)
,f_lim_max(f_lim_max)
,gamma_min(gamma_min)
,gamma_max(gamma_max)
,x_min(Data::get_instance().get_x_min())
,x_max(Data::get_instance().get_x_max())
,y_min(Data::get_instance().get_y_min())
,y_max(Data::get_instance().get_y_max())
{

}

void MyConditionalPrior::from_prior(DNest4::RNG& rng)
{
  f_lim = exp(log(f_lim_min) + log(f_lim_max/f_lim_min)*rng.rand());
  gamma = gamma_min + (gamma_max-gamma_min)*rng.rand();
}

double MyConditionalPrior::perturb_hyperparameters(DNest4::RNG& rng)
{
  double logH = 0.;

  int which = rng.rand_int(2);

  if(which == 0)
  {
    // perturb flux limit
    double log_f_lim = log(f_lim);
    log_f_lim += log(f_lim_max/f_lim_min)*rng.randh();
    DNest4::wrap(log_f_lim, log(f_lim_min), log(f_lim_max));
    f_lim = exp(log_f_lim);
    // logH += 0;
  }
  else if(which == 1)
  {
    // perturb flux slope
    gamma += (gamma_max-gamma_min)*rng.randh();
    DNest4::wrap(gamma, gamma_min, gamma_max);
    // logH += 0;
  }

  return logH;
}

double MyConditionalPrior::log_pdf(const vector<double>& vec) const
{
  Uniform Px(x_min, x_max);
  Uniform Py(y_min, y_max);
  Pareto Pflux(f_lim, gamma);

  double logp = 0.;

  logp += Px.log_pdf(vec[0]);
  logp += Py.log_pdf(vec[1]);
  logp += Pflux.log_pdf(vec[2]);
  
  return logp;
}

void MyConditionalPrior::from_uniform(vector<double>& vec) const
{
  Uniform Px(x_min, x_max);
  Uniform Py(y_min, y_max);
  Pareto Pflux(f_lim, gamma);

  vec[0] = Px.cdf_inverse(vec[0]);
  vec[1] = Py.cdf_inverse(vec[1]);
  vec[2] = Pflux.cdf_inverse(vec[2]); 
}

void MyConditionalPrior::to_uniform(vector<double>& vec) const
{
  Uniform Px(x_min, x_max);
  Uniform Py(y_min, y_max);
  Pareto Pflux(f_lim, gamma);

  vec[0] = Px.cdf(vec[0]);
  vec[1] = Py.cdf(vec[1]);
  vec[2] = Pflux.cdf(vec[2]); 
}

void MyConditionalPrior::print(ostream& out) const
{
  out << f_lim << " " << gamma << " ";
}
