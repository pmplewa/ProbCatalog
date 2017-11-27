#ifndef DNest4_RJObject_TestCase_PSF
#define DNest4_RJObject_TestCase_PSF

#include <vector>
#include <ostream>
#include "../DNest4/code/DNest4.h"

class GaussianPSF
{
  private:
    double scale; // image scale (PSF size)
    double scale_min, scale_max;

  public:
    GaussianPSF(double scale_min, double scale_max);

    double evaluate(double x, double y) const;

    void from_prior(DNest4::RNG& rng);
    double perturb(DNest4::RNG& rng);

    double log_pdf() const { return 0; };

    void print(std::ostream& out) const;
};

class DoubleGaussianPSF
{
  private:
    double scale; // size of PSF core
    double scale_min, scale_max;

    double scale2; // size of PSF halo
    double scale2_min, scale2_max;

    double frac; // core/halo fraction
    double frac_min, frac_max;

  public:
    DoubleGaussianPSF(
        double scale_min, double scale_max,
        double scale2_min, double scale2_max,
        double frac_min, double frac_max);

    double evaluate(double x, double y) const;

    void from_prior(DNest4::RNG& rng);
    double perturb(DNest4::RNG& rng);

    double log_pdf() const;

    void print(std::ostream& out) const;
};

class MoGPSF
{
  private:
    int n_comp; // number of components
    int n_param; // number of parameters

    double scale; // image scale parameter
    double scale_min, scale_max;

    std::vector<double> frac; // amplitudes
    std::vector<std::vector<double>> mean; // mean vectors
    std::vector<std::vector<double>> cov; // covariance matrices

  public:
    MoGPSF(
      int n_comp, const char* param_file,
      double scale_min, double scale_max);

    double evaluate(double x, double y) const;

    void from_prior(DNest4::RNG& rng);
    double perturb(DNest4::RNG& rng);

    double log_pdf() const { return 0; };

    void print(std::ostream& out) const;
};

class ShapeletPSF
{
  private:
    int n_max; // maximum order
    int n_param; // number of parameters (nx + ny < n_max)

    double scale; // image scale (PSF size)
    double scale_min, scale_max;

    // pre-computed coefficients
    std::vector<double> coeff1;
    std::vector<double> coeff2;

    std::vector<double> param; // parameter vector

    // parameter indexing
    int index(int nx, int ny) const;

    // basis functions
    double phi1d(double x, int n) const;
    //double B1d(double x, int n) const;
    double phi2d(double x, double y, int nx, int ny) const;
    double B2d(double x, double y, int nx, int ny) const;

  public:
    ShapeletPSF(
      int n_max, const char* param_file,
      double scale_min, double scale_max);

    double evaluate(double x, double y) const;

    void from_prior(DNest4::RNG& rng);
    double perturb(DNest4::RNG& rng);

    double log_pdf() const { return 0; };

    double get_flux();
    void normalize();

    void print(std::ostream& out) const;
};

#endif
