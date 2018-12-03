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
  int model_;
  double const lambda = -0.1;

  double f(double t,const double &y) const
  {
    if(model_ == 1)
    {
    return 1-cos(t)*(cos(t)-1) - y*y;
    }
    else if (model_ == 2)
    {
      return 1-cos(t)*(cos(t)-exp(lambda*t)) - exp(-2*lambda*t)*y*y + lambda*y;;
    }
    else
    {
      return -y;
    }
  }
  double df(double t,const double &y) const
  {
    if(model_ == 1)
    {
    return -2*y;
    }
    else if (model_ == 2)
    {
      return -2*exp(-2*lambda*t)*y + lambda;
    }
    else
    {
      return -1;
    }
  }
  double T() const
  {
    if (model_ == 1)
    {
    return 2*M_PI;
    }
    else if (model_ == 2)
    {
      return 10;
    }
    else
    {
      return 1;
    }
  }
  double y0() const
  {
    if(model_ == 1)
    {
    return 0.;
    }
    else if (model_ == 2)
    {
      return 0.;
    }
    else
    {
      return 1.;
    }
  }
  double exact(double t) const
  {
    if(model_ == 1)
    {
    return sin(t);
    }
    else if(model_ == 2)
    {
      return exp(lambda*t)*sin(t);
    }
    else
    {
      return exp(-t);
    }
  }
  };


#endif // MODELS_HH
