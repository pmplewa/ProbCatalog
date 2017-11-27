#include <iostream>
#include "DNest4/code/DNest4.h"
#include "MyModel.h"
#include "Data.h"

using namespace std;
using namespace DNest4;

int main(int argc, char** argv)
{
  static_assert(numeric_limits<float>::is_iec559, "ERROR: IEEE 754 required");

  Data::get_instance().load(
    "data/test_metadata.txt",
    "data/test_image.txt");

  Sampler<MyModel> sampler = setup<MyModel>(argc, argv);
  sampler.run();

  return 0;
}
