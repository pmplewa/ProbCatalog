#ifndef DNest4_TestCase_Lookup
#define DNest4_TestCase_Lookup

#include <vector>

class Lookup
{
  private:
    int n_step;
    double x_min, x_max;

    double dx;
    
    int n_max;

    std::vector<std::vector<double>> _H;

    Lookup();

    static Lookup instance;

  public:
    static double evaluate(double x, int n);
};

#endif
