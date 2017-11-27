#ifndef DNest4_TestCase_MyModel
#define DNest4_TestCase_MyModel

#include <vector>
#include "DNest4/code/DNest4.h"
#include "MyConditionalPrior.h"
#include "PSF.h"
#include "Noise.h"

class MyModel
{
  private:
    // model image
    int ni, nj; // image dimensions
    std::vector<std::vector<long double>> image;
    void calculate_image();

    // point sources in the image
    DNest4::RJObject<MyConditionalPrior> objects;

    // PSF model
    #ifdef PSFModel_GaussianPSF
    GaussianPSF psf; 
    #elif PSFModel_DoubleGaussianPSF
    DoubleGaussianPSF psf;
    #elif PSFModel_MoGPSF
    MoGPSF psf;
    #elif PSFModel_ShapeletPSF
    ShapeletPSF psf;    
    #else
    #error "# ERROR: Invalid PSF model"
    #endif    

    // noise model
    #ifdef NoiseModel_UniformNoise 
    UniformNoise noise;
    #elif NoiseModel_VariableNoise
    VariableNoise noise;
    #else
    #error "# ERROR: Invalid noise model"
    #endif

  public:
    MyModel();

    // generate particles from the prior
    void from_prior(DNest4::RNG& rng);

    // propose Metropolis-Hastings steps
    double perturb(DNest4::RNG& rng);

    // calculate the likelihood
    double log_likelihood();

    // print output
    void print(std::ostream& out) const;
    std::string description() const;
};

#endif
