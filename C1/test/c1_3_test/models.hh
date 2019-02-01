#ifndef MODELS_HH
#define MODELS_HH

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>

struct Test
{
  double f(double t,const double &y) const
  {
    return -y;
  }
  double df(double t,const double &y) const
  {
    return -1;
  }
  double T() const 
  {
    return 1.;
  }
  double y0() const 
  {
    return 1.;
  }
  double exact(double t) const 
  {
    return exp(-t);
  }
};

#endif // MODELS_HH
