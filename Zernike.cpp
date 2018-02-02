#include "Zernike.h"
#include "DNest4/code/DNest4.h"
#include <math.h>
#include <stdexcept>
#include <utility>
#include <boost/math/special_functions/factorials.hpp>

using namespace std;
using namespace boost::math;

Zernike::Zernike(int j_max, double param_min, double param_max)
:j_max(j_max)
,param(j_max+1, 0.)
,param_min(param_min)
,param_max(param_max)
{
  int n, m, abs_m;
  int n_terms;

  for(int j=0; j<=j_max; j++)
  {
    tie(n, m) = index(j);
    abs_m = abs(m);

    n_terms = (n - abs_m)/2 + 1;
    R_norm.push_back(m == 0 ? sqrt(n + 1) : sqrt(2*(n + 1)));
    R_order.push_back(vector<double>(n_terms));
    R_coeff.push_back(vector<double>(n_terms));
    for(int k=0; k<n_terms; k++)
    {
      R_order[j][k] = n - 2*k;
      if((n - abs_m) % 2 == 0)
        R_coeff[j][k] = pow(-1, k)
                      * factorial<double>(n - k)
                      / factorial<double>(k)
                      / factorial<double>((n + abs_m)/2 - k)
                      / factorial<double>((n - abs_m)/2 - k);
      else
        R_coeff[j][k] = 0.;
    }
  }
}

std::pair<int, int> Zernike::index(int j) const
{
  if(j < 0)
    throw invalid_argument("# ERROR: allowed parameters must fulfill j >= 0");

  int n, m;

  n = 0;
  while(j > n)
  {
    n += 1;
    j -= n;
  }
  m = 2*j - n;

  return make_pair(n, m);
}

double Zernike::evaluate(double x, double y) const
{
  double value = 0.;

  double r = sqrt(pow(x, 2) + pow(y, 2));
  double phi = atan2(y, x);

  int n, m, abs_m;
  int n_terms;

  for(int j=0; j<=j_max; j++)
  {
    tie(n, m) = index(j);
    abs_m = abs(m);

    n_terms = (n - abs_m)/2 + 1;
    for(int k=0; k<n_terms; k++)
      if(m < 0)
        value += param[j] * R_norm[j] * R_coeff[j][k]
               * pow(r, R_order[j][k]) * sin(abs_m * phi);
      else
        value += param[j] * R_norm[j] * R_coeff[j][k]
               * pow(r, R_order[j][k]) * cos(abs_m * phi);        
  }

  return value;
}

void Zernike::from_prior(DNest4::RNG& rng)
{
  for(int j=0; j<=j_max; j++)
    param[j] = param_min + (param_max - param_min)*rng.rand();
}

double Zernike::perturb(DNest4::RNG& rng)
{
  double logH = 0.;

  int which = rng.rand_int(j_max+1);

  param[which] += (param_max - param_min)*rng.randh();
  DNest4::wrap(param[which], param_min, param_max);
  // logH += 0;

  return logH;
}

void Zernike::print(std::ostream& out) const
{
  for(int j=0; j<=j_max; j++)
    out << param[j] << " ";
}
