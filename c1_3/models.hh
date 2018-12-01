#ifndef MODELS_HH
#define MODELS_HH

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <math.h>

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

struct Test1
{
  double f(double t,const double &y) const
  {
    return 1-cos(t)*(cos(t)-1) - y*y;
  }
  double df(double t,const double &y) const
  {
    return -2*y;
  }
  double T() const
  {
    return 2*M_PI;
  }
  double y0() const
  {
    return 0.;
  }
  double exact(double t) const
  {
    return sin(t);
  }
};

#endif // MODELS_HH
