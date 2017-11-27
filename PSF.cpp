#include "PSF.h"
#include "DNest4/code/DNest4.h"
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <cmath>
#include <boost/math/special_functions/binomial.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/hermite.hpp>
//#include "Lookup.h"

using namespace std;
using namespace boost::math;

GaussianPSF::GaussianPSF(double scale_min, double scale_max)
:scale_min(scale_min)
,scale_max(scale_max)
{

}

double GaussianPSF::evaluate(double x, double y) const
{
  double r_2 = pow(x, 2) + pow(y, 2);
  double scale_2 = pow(scale, 2);

  return exp(-0.5 * r_2/scale_2)/(2.*M_PI * scale_2);
}

void GaussianPSF::from_prior(DNest4::RNG& rng)
{
  scale = exp(log(scale_min) + log(scale_max/scale_min)*rng.rand());
}

double GaussianPSF::perturb(DNest4::RNG& rng)
{
  double logH = 0.;

  // perturb image scale
  double log_scale = log(scale);
  log_scale += log(scale_max/scale_min)*rng.randh();
  DNest4::wrap(log_scale, log(scale_min), log(scale_max));
  scale = exp(log_scale);
  // logH += 0;

  return logH;
}

void GaussianPSF::print(std::ostream& out) const
{
  out << scale << " ";
}

DoubleGaussianPSF::DoubleGaussianPSF(
  double scale_min, double scale_max,
  double scale2_min, double scale2_max,
  double frac_min, double frac_max)
:scale_min(scale_min)
,scale_max(scale_max)
,scale2_min(scale2_min)
,scale2_max(scale2_max)
,frac_min(frac_min)
,frac_max(frac_max)
{

}

double DoubleGaussianPSF::evaluate(double x, double y) const
{
  double value = 0;

  double r_2 = pow(x, 2) + pow(y, 2);
  double scale_2 = pow(scale, 2);
  double scale2_2 = pow(scale2, 2);

  value += frac * exp(-0.5 * r_2/scale_2)/(2.*M_PI * scale_2);
  value += (1. - frac) * exp(-0.5 * r_2/scale2_2)/(2.*M_PI * scale2_2);

  return value;
}

void DoubleGaussianPSF::from_prior(DNest4::RNG& rng)
{
  scale = exp(log(scale_min) + log(scale_max/scale_min)*rng.rand());
  scale2 = exp(log(scale2_min) + log(scale2_max/scale2_min)*rng.rand());
  frac = frac_min + (frac_max - frac_min)*rng.rand();
}

double DoubleGaussianPSF::perturb(DNest4::RNG& rng)
{
  double logH = 0.;

  int which = rng.rand_int(3);

  if(which == 0)
  {
    // perturb size of PSF core
    double log_scale = log(scale);
    log_scale += log(scale_max/scale_min)*rng.randh();
    DNest4::wrap(log_scale, log(scale_min), log(scale_max));
    scale = exp(log_scale);
    // logH += 0;
  }
  else if(which == 1)
  {
    // perturb size of PSF halo
    double log_scale2 = log(scale2);
    log_scale2 += log(scale2_max/scale2_min)*rng.randh();
    DNest4::wrap(log_scale2, log(scale2_min), log(scale2_max));
    scale2 = exp(log_scale2);
    // logH += 0;
  }   
  else if(which == 2)
  {
    // perturb core/halo fraction
    frac += (frac_max - frac_min)*rng.randh();
    DNest4::wrap(frac, frac_min, frac_max);
    // logH += 0;  
  }     

  return logH;
}

double DoubleGaussianPSF::log_pdf() const
{
  if(scale2 < scale)
    return -numeric_limits<double>::max();

  return 0.;
}

void DoubleGaussianPSF::print(std::ostream& out) const
{
  out << scale << " " << scale2 << " " << frac << " ";
}

MoGPSF::MoGPSF(
  int n_comp, const char* param_file,
  double scale_min, double scale_max)
:n_comp(n_comp)
,n_param(6*n_comp)
,scale_min(scale_min)
,scale_max(scale_max)
,frac(n_param)
,mean(n_param, vector<double>(2))
,cov(n_param, vector<double>(3))
{
  // load parameter values
  fstream fin(param_file, ios::in);
  if(!fin)
  {
    cerr << "# ERROR: could not open file " << param_file << endl;
    abort();
  }
  for(int n=0; n<n_comp; n++)
    fin >> frac[n];
  for(int n=0; n<n_comp; n++)
    for(int i=0; i<2; i++)
      fin >> mean[n][i];
  for(int n=0; n<n_comp; n++)
    for(int i=0; i<3; i++)
      fin >> cov[n][i];
  fin.close();
}

double MoGPSF::evaluate(double x, double y) const
{
  double value = 0;

  double x0, y0;
  double sx_2, sy_2, sxy_2;
  double scale_2 = pow(scale, 2);

  for(int n=0; n<n_comp; n++)
  {
    x0 = scale * mean[n][0];
    y0 = scale * mean[n][1];

    sx_2 = scale_2 * cov[n][0];
    sy_2 = scale_2 * cov[n][1];
    sxy_2 = scale_2 * cov[n][2];

    value += frac[n]
        * exp(-0.5*(
              pow(x-x0, 2)/sx_2
            - 2.*sxy_2*(x-x0)*(y-y0)/(sx_2*sy_2)
            + pow(y-y0, 2)/sy_2
          )/(1.-pow(sxy_2, 2)/(sx_2*sy_2)))
        /(2.*M_PI * sqrt(sx_2*sy_2-pow(sxy_2, 2)));
  }

  return value;
}

void MoGPSF::from_prior(DNest4::RNG& rng)
{
  scale = exp(log(scale_min) + log(scale_max/scale_min)*rng.rand());
}

double MoGPSF::perturb(DNest4::RNG& rng)
{
  double logH = 0.;

  // perturb image scale
  double log_scale = log(scale);
  log_scale += log(scale_max/scale_min)*rng.randh();
  DNest4::wrap(log_scale, log(scale_min), log(scale_max));
  scale = exp(log_scale);
  // logH += 0;

  return logH;
}

void MoGPSF::print(std::ostream& out) const
{
  out << scale << " ";
}

ShapeletPSF::ShapeletPSF(
  int n_max, const char* param_file,
  double scale_min, double scale_max)
:n_max(n_max)
,n_param(n_max+(n_max*n_max-n_max)/2)
,scale_min(scale_min)
,scale_max(scale_max)
,coeff1(n_max)
,coeff2(n_param)
,param(n_param)
{
  // pre-compute coefficients for shapelet evaluation
  for(int n=0; n<n_max; n++)
    coeff1[n] = pow(pow(2, n) * pow(M_PI, 0.5) * factorial<double>(n), -0.5);

  // pre-compute coefficients for flux calculation
  for(int nx=0; nx<n_max; nx+=2)
    for(int ny=0; ny<n_max; ny+=2)
      if(nx+ny < n_max)
        coeff2[index(nx, ny)] = pow(M_PI, 0.5) * pow(2., 0.5*(2-nx-ny))
          * pow(binomial_coefficient<double>(nx, nx/2), 0.5)
          * pow(binomial_coefficient<double>(ny, ny/2), 0.5);

  // set up parameter priors
  fstream fin(param_file, ios::in);
  if(!fin)
  {
    cerr << "# ERROR: could not open file " << param_file << endl;
    abort();
  }
  for(int n=0; n<n_param; n++)
    fin >> param[n];
  fin.close();
}

int ShapeletPSF::index(int nx, int ny) const
{
  //if(!(nx+ny < n_max))
  //  throw invalid_argument("# ERROR: maximum order exceeded");
  return n_max - ny - (nx*(1 + nx - 2*n_max))/2 - 1;
}

double ShapeletPSF::phi1d(double x, int n) const
{
  return coeff1[n] * hermite(n, x) * exp(-pow(x, 2)/2);
  //return coeff1[n] * Lookup::evaluate(x, n);
}

//double B1d(double x, int n, double scale) const
//{
//  return phi1d(x/scale, n)/sqrt(scale)
//}

double ShapeletPSF::phi2d(double x, double y, int nx, int ny) const
{
  return phi1d(x, nx) * phi1d(y, ny);
}

double ShapeletPSF::B2d(double x, double y, int nx, int ny) const
{
  return phi2d(x/scale, y/scale, nx, ny)/scale;
}

double ShapeletPSF::evaluate(double x, double y) const
{
  double value = 0;

  for(int nx=0; nx<n_max; nx++)
    for(int ny=0; ny<n_max; ny++)
      if(nx+ny < n_max)
        value += param[index(nx, ny)] * B2d(x, y, nx, ny);

  return value;
}

void ShapeletPSF::from_prior(DNest4::RNG& rng)
{
  scale = exp(log(scale_min) + log(scale_max/scale_min)*rng.rand());

  normalize();
}

double ShapeletPSF::perturb(DNest4::RNG& rng)
{
  double logH = 0.;

  // perturb image scale
  double log_scale = log(scale);
  log_scale += log(scale_max/scale_min)*rng.randh();
  DNest4::wrap(log_scale, log(scale_min), log(scale_max));
  scale = exp(log_scale);
  // logH += 0;

  normalize();

  return logH;
}

double ShapeletPSF::get_flux()
{
  double flux = 0.;

  for(int nx=0; nx<n_max; nx+=2)
    for(int ny=0; ny<n_max; ny+=2)
      if(nx+ny < n_max)
        flux += coeff2[index(nx, ny)] * param[index(nx, ny)];

  flux *= scale;

  return flux;
}

void ShapeletPSF::normalize()
{
  double flux = get_flux();
  for(int n=0; n<n_param; n++)
    param[n] /= flux;
}

void ShapeletPSF::print(std::ostream& out) const
{
  out << scale << " ";
  for(int n=0; n<n_param; n++)
    out << param[n] << " ";
}

