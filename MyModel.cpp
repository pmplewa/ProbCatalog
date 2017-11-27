#include "MyModel.h"
#include "DNest4/code/DNest4.h"
#include "Data.h"
#include <cmath>

using namespace std;

MyModel::MyModel()
:ni(Data::get_instance().get_ni())
,nj(Data::get_instance().get_nj())
,image(ni, vector<long double>(nj))
,objects(
  3, // num_dimensions
  100, // max_num_components
  false, // fixed
  MyConditionalPrior(
    1E-6, // f_lim_min
    1E+0, // f_lim_max
    0., // gamma_min
    100. // gamma_max
  ),
  DNest4::PriorType::uniform
)
,psf(
  #ifdef PSFModel_GaussianPSF
  GaussianPSF(
    1E-3, // scale_min
    1E+0 // scale_max
  )
  #elif PSFModel_DoubleGaussianPSF
  DoubleGaussianPSF(
    1E-3, // scale_min
    1E+0, // scale_max
    1E-3, // scale2_min
    1E+0, // scale2_max  
    0., // frac_min
    1. // frac_max
  )
  #elif PSFModel_MoGPSF
  MoGPSF(
    20, // n_comp
    "data/test_param.txt", // param_file
    1E-3, // scale_min
    1E+0 // scale_max
  )
  #elif PSFModel_ShapeletPSF
  ShapeletPSF(
    15, // n_max
    "data/test_param.txt", // param_file
    1E-3, // scale_min
    1E+0 // scale_max
  )  
  #endif
)
,noise(
  #ifdef NoiseModel_UniformNoise 
  UniformNoise(
    -1., // mean_min
    1., // mean_max
    1E-6, // sigma_min
    1E+0 // sigma_max
  )
  #elif NoiseModel_VariableNoise
  VariableNoise(
    6, // j_max
    -1., // param_min
    1., // param_max
    1E-6, // sigma_min
    1E+0 // sigma_max    
  )
  #endif
)
{

}

void MyModel::from_prior(DNest4::RNG& rng)
{
  objects.from_prior(rng);
  noise.from_prior(rng);
  psf.from_prior(rng);
}

double MyModel::perturb(DNest4::RNG& rng)
{
  double logH = 0.;

  int which = rng.rand_int(2);

  if(which == 0)
  {
    // perturb objects
    logH += objects.perturb(rng);
  }
  else if(which == 1)
  {
    which = rng.rand_int(2);
    
    if (which == 0)
    {
      // perturb noise (including background)
      logH += noise.perturb(rng);
    }
    else if(which == 1)
    {
      // perturb PSF (including image scale)
      logH += psf.perturb(rng);
    }
  }

  return logH;
}

void MyModel::calculate_image()
{
  // get objects
  const vector<vector<double>>& components = objects.get_components();

  // get image coordinates
  const vector<vector<double>>& x = Data::get_instance().get_x_rays();
  const vector<vector<double>>& y = Data::get_instance().get_y_rays();

  // reset the image
  for(int i=0; i<ni; i++)
    for(int j=0; j<nj; j++)
        image[i][j] = noise.get_mean(x[i][j], y[i][j]);     

  double xc, yc, f;      

  // plant sources in the image
  for(size_t k=0; k<components.size(); k++)
  {
    xc = components[k][0];
    yc = components[k][1];
    f = components[k][2];

    for(int i=0; i<ni; i++)
      for(int j=0; j<nj; j++)
        image[i][j] += f * psf.evaluate(x[i][j] - xc, y[i][j] - yc);
  }
}

double MyModel::log_likelihood()
{
  // get data
  const vector<vector<double>>& data = Data::get_instance().get_image();
  const vector<vector<double>>& mask = Data::get_instance().get_mask();

  double logL = 0.;

  double sigma = noise.get_sigma();
  double sigma_2 = sigma*sigma;

  calculate_image();

  for(int i=0; i<ni; i++)
    for(int j=0; j<nj; j++)
      if(mask[i][j])
        logL += -0.5*log(2.*M_PI*sigma_2)
          -0.5*pow(data[i][j] - image[i][j], 2)/sigma_2;

  logL += psf.log_pdf();

  return logL;
}

void MyModel::print(std::ostream& out) const
{
  out << setprecision(6);
  for(int i=0; i<ni; i++)
    for(int j=0; j<nj; j++)
      out << image[i][j] << " ";
  out << setprecision(10);
  objects.print(out);
  noise.print(out);
  psf.print(out);
}

string MyModel::description() const
{
  return string("objects");
}
