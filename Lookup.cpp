#include "Lookup.h"
#include <cmath>
#include <boost/math/special_functions/hermite.hpp>

using namespace std;
using namespace boost::math;

Lookup Lookup::instance;

Lookup::Lookup()
:n_step(20001)
,x_min(-10.)
,x_max(10.)
,dx((x_max-x_min)/(n_step-1))
,n_max(15)
{
  _H.assign(n_max, vector<double>(n_step));

  double x;

  for(int n=0; n<n_max; n++)
    for(int i=0; i<n_step; i++)
    {
      x = x_min + i*dx;
      _H[n][i] = hermite(n, x) * exp(-pow(x, 2)/2);
    }
}

double Lookup::evaluate(double x, int n)
{
  int i = static_cast<int>((x-instance.x_min)/instance.dx);
  double frac = (x - (instance.x_min + i*instance.dx))/instance.dx;

  if(i < 0 || i >= instance.n_step || n >= instance.n_max)
    return 0.;

  return frac*instance._H[n][i+1] + (1. - frac)*instance._H[n][i];  
}
