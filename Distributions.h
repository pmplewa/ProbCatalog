#ifndef DNest4_TestCase_Distributions
#define DNest4_TestCase_Distributions

#include "../DNest4/code/Distributions/ContinuousDistribution.h"

bool is_close(double x, double y, double eps=1e-6);

// Uniform Distribution
class Uniform:public DNest4::ContinuousDistribution
{
  private:
    double min, max;

  public:
    Uniform(double min, double max);

    double log_pdf(double x) const;

    double cdf(double x) const;
    double cdf_inverse(double x) const;
};

// Pareto Distribution
class Pareto:public DNest4::ContinuousDistribution
{
  private:
    double min, gamma;

  public:
    Pareto(double min, double gamma);

    double log_pdf(double x) const;

    double cdf(double x) const;
    double cdf_inverse(double x) const;
};

// Truncated Pareto Distribution
class ParetoTrunc:public DNest4::ContinuousDistribution
{
  private:
    double min, max, gamma;

  public:
    ParetoTrunc(double min, double max, double gamma);

    double log_pdf(double x) const;
    
    double cdf(double x) const;
    double cdf_inverse(double x) const;
};

#endif
