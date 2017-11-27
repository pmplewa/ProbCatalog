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
,param(j_max, 0.)
,param_min(param_min)
,param_max(param_max)
{
  int n, m, j;
  int n_terms;

  for(int i=0; i<j_max; i++)
  {
    j = i + 1;
    tie(n, m) = idx(j);
    m = abs(m);

    n_terms = (n - m)/2 + 1;
    R_norm.push_back(m == 0 ? sqrt(n + 1) : sqrt(2*(n + 1)));
    R_order.push_back(vector<double>(n_terms));
    R_coeff.push_back(vector<double>(n_terms));
    for(int k=0; k<n_terms; k++)
    {
      R_order[i][k] = n - 2*k;
      if((n - m) % 2 == 0)
        R_coeff[i][k] = pow(-1, k)
                      * factorial<double>(n - k)
                      / factorial<double>(k)
                      / factorial<double>((n + m)/2 - k)
                      / factorial<double>((n - m)/2 - k);
      else
        R_coeff[i][k] = 0.;
    }
  }
}

std::pair<int, int> Zernike::idx(int j) const
{
  if(j < 1)
    throw invalid_argument("# ERROR: allowed parameters must fulfill j >= 1");

  int n, m, r;

  n = 0;
  r = j - 1;
  while(r > n)
  {
    n += 1;
    r -= n;
  }

  m = pow(-1, j) * ((n % 2) + 2*int((r + ((n + 1) % 2)) / 2.));

  return make_pair(n, m);
}

double Zernike::evaluate(double x, double y) const
{
  double value = 0.;

  double r = sqrt(pow(x, 2) + pow(y, 2));
  double phi = atan2(y, x);

  int n, m, j;
  int n_terms;

  for(int i=0; i<j_max; i++)
  {
    j = i + 1;
    tie(n, m) = idx(j);
    m = abs(m);

    n_terms = (n - m)/2 + 1;
    for(int k=0; k<n_terms; k++)
      if((j % 2) == 0)
        value += param[i] * R_norm[i] * R_coeff[i][k]
               * pow(r, R_order[i][k]) * cos(m * phi);
      else
        value += param[i] * R_norm[i] * R_coeff[i][k]
               * pow(r, R_order[i][k]) * sin(m * phi);
  }

  return value;
}

void Zernike::from_prior(DNest4::RNG& rng)
{
  for(int i=0; i<j_max; i++)
    param[i] = param_min + (param_max - param_min)*rng.rand();
}

double Zernike::perturb(DNest4::RNG& rng)
{
  double logH = 0.;

  int which = rng.rand_int(j_max);

  param[which] += (param_max - param_min)*rng.randh();
  DNest4::wrap(param[which], param_min, param_max);
  // logH += 0;

  return logH;
}

void Zernike::print(std::ostream& out) const
{
  for(int i=0; i<j_max; i++)
    out << param[i] << " ";
}
