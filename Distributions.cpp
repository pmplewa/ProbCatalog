#include "Distributions.h"
#include <stdexcept>
#include <iostream>

using namespace std;

bool is_close(double x, double y, double eps)
{
  return abs(x - y) <= eps * abs(x);
}

Uniform::Uniform(double min, double max)
:min(min)
,max(max)
{
  if(max <= min)
    throw domain_error("# ERROR: allowed parameters must fulfill max > min");
}

double Uniform::cdf(double x) const
{
  if(x < min)
    return 0.;
  else if(x > max)
    return 1.;
  else
    return (x - min)/(max - min);
}

double Uniform::cdf_inverse(double x) const
{
  if(x < 0. || x > 1.)
      throw domain_error("# ERROR: input to cdf_inverse must be in [0, 1]");

  return min + (max - min)*x;
}

double Uniform::log_pdf(double x) const
{
  if(x < min)
    return -numeric_limits<double>::infinity();
  else if(x > max)
    return -numeric_limits<double>::infinity();
  else
    return log(1./(max - min));
}

Pareto::Pareto(double min, double gamma)
:min(min)
,gamma(gamma)
{
  if(min <= 0.)
    throw domain_error("# ERROR: allowed parameters must fulfill min > 0");

  if(gamma <= 0.)
    throw domain_error("# ERROR: allowed parameters must fulfill gamma > 0");
}

double Pareto::cdf(double x) const
{
  double alpha = 1./gamma;

  if(x < min)
    return 0.;
  else
    return 1. - pow(min/x, alpha);
}

double Pareto::cdf_inverse(double x) const
{
  if(x < 0. || x > 1.)
      throw domain_error("# ERROR: input to cdf_inverse must be in [0, 1]");

  return min * pow(1. - x, -gamma); // cdf_inverse(1) = inf
}

double Pareto::log_pdf(double x) const
{
  double alpha = 1./gamma;

  if(x < min)
    return -numeric_limits<double>::infinity();
  else
    return log(alpha) + alpha*log(min) - (alpha + 1.)*log(x);
}

ParetoTrunc::ParetoTrunc(double min, double max, double gamma)
:min(min)
,max(max)
,gamma(gamma)
{
  if(min <= 0.)
    throw domain_error("# ERROR: allowed parameters must fulfill min > 0");

  if(max <= min)
    throw domain_error("# ERROR: allowed parameters must fulfill max > min");

  if(gamma <= 0.)
    throw domain_error("# ERROR: allowed parameters must fulfill gamma > 0");
}

double ParetoTrunc::cdf(double x) const
{
  double alpha = 1./gamma;

  if(x <= min)
    return 0.;
  else if(x >= max)
    return 1.;
  else
    return (pow(min/x, alpha) - 1.)/(pow(min/max, alpha) - 1.);
}

double ParetoTrunc::cdf_inverse(double x) const
{
  if(x < 0. || x > 1.)
      throw domain_error("# ERROR: input to cdf_inverse must be in [0, 1]");

  double alpha = 1./gamma;

  if(is_close(x, 0.))
    return min;
  else if(is_close(x, 1.))
    return max;
  else
    return pow(pow(min, -alpha) * (x*(pow(min/max, alpha) - 1.) + 1.), -gamma);
}

double ParetoTrunc::log_pdf(double x) const
{
  double alpha = 1./gamma;

  if(x <= min)
    return -numeric_limits<double>::infinity();
  else if(x >= max)
    return -numeric_limits<double>::infinity();
  else
    return log(alpha) + alpha*log(min) - (alpha + 1.)*log(x) - log(1. - pow(min/max, alpha));
}
