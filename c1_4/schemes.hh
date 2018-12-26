#ifndef SCHEMES_HH
#define SCHEMES_HH

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>

class DIRK
{
public:

  DIRK(int stages)
  : stages_(stages),
    a_(stages*stages), b_(stages), c_(stages)
  {
    for (int i=0;i<stages_;++i)
      for (int j=0;j<stages_;++j)
        a_[i*stages_+j] = 0.;
     // putting in forward Euler here (although we are not using them at the moment)
     b_[0] = 1.;
     c_[0] = 0.;
  }

  // evolve the solution y(t) to y(t+h)
  // take y,t,h, and the model and return the approximation at the next
  // time level
  // The Model is given by a template argument, which has to provide the
  // following methods:
  // - model.f(t,y)
  // - model.df(t,y)
  // - model.N
    template <class Model>
  Vector evolve(const Vector &y, double t, double h, const Model &model) const
  {
    Vector ret = y;
    std::vector<Vector> k(stages_, Vector(model.N_));

	for (int s=0; s<stages_; ++s )
  {
	  // compute k_s and store inside k[s]
    Vector sum = y;
    for(int j = 0; j < s; j++)
    {
      sum = sum + k[j]*(h*a(s,j));
    }
	  if (a(s,s) == 0) // explicit case
    {
      if(c_[s] != 0)
      {
        sum = model.f(c_[s]*h,sum);
        k[s] = model.f(1,sum);
      }
      else
      {
        k[s] = model.f(1,sum);
      }
    }
    else
    {
		  // implicit case
		  // solve via Newton method
      std::cout << "We'll deal with this later" << std::endl;
	  }
	  // Increment the return value by the current k[s]
	  ret =ret + k[s]*(h*b_[s]);
	}
    return ret;
  }

protected:

  const double a(int i, int j) const { return a_[i*stages_+j]; }
  double& a(int i, int j) { return a_[i*stages_+j]; }

  int stages_;

  std::vector<double> a_,b_,c_;
};


class FE: public DIRK
{
public:

  FE() : DIRK(1)
  {
  	a(0,0) = 0.;
  	b_[0] = 1.;
  	c_[0] = 0.;
  }
};

#endif

class Heun3: public DIRK
{
public:

  Heun3() : DIRK(3)
  {
	a(1,0) = double(1)/double(3);
  a(2,1) = double(2)/double(3);
	b_[0] = 0.25;
  b_[1] = 0.;
  b_[2] = 0.75;
  c_[0] = 0.;
  c_[1] = double(1)/double(3);
  c_[2] = double(2)/double(3);
  }
};
