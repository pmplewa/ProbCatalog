#include "Data.h"
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

Data Data::instance;

Data::Data()
:x_min(-1.)
,x_max(1.)
,y_min(-1.)
,y_max(1.)
,r_max(1.)
{

}

void Data::load(const char* metadata_file, const char* image_file)
{
  // load the metadata
  fstream fin(metadata_file, ios::in);
  if(!fin)
  {
    cerr << "# ERROR: could not open file " << metadata_file << endl;
    abort();
  }
  fin >> ni >> nj;
  fin.close();

  // calculate pixel widths
  dx = (x_max - x_min)/ni;
  dy = (y_max - y_min)/nj;

  // check that pixels are square
  if(abs(log(dx/dy)) >= 1e-6)
  {
    cerr << "# ERROR: pixels are not square" << endl;
    abort();
  }

  // calculate image coordinates
  x_rays.assign(ni, vector<double>(nj));
  y_rays.assign(ni, vector<double>(nj));
  for(int i=0; i<ni; i++)
    for(int j=0; j<nj; j++)
    {
      x_rays[i][j] = x_min + (j + 0.5)*dx;
      y_rays[i][j] = y_max - (i + 0.5)*dy;
    }

  // load the image data
  fin.open(image_file, ios::in);
  if(!fin)
  {
    cerr << "# ERROR: could not open file " << image_file << endl;
    abort();
  }
  image.assign(ni, vector<double>(nj));
  mask.assign(ni, vector<double>(nj));
  for(int i=0; i<ni; i++)
    for(int j=0; j<nj; j++)
    {
      fin >> image[i][j];
      if(sqrt(pow(x_rays[i][j], 2) + pow(y_rays[i][j], 2)) < r_max)
        mask[i][j] = 1.; // inside
      else
        mask[i][j] = 0.; // outside
    }
  fin.close();
}
