#ifndef SCHEMES_HH
#define SCHEMES_HH

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>
#include "models.hh"


class DIRK
{
public:

  DIRK(int stages) // Default constructor
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
  // takes y,t,h, and the model and returns the approximation at the next
  // time level
  // The Model is given by a template argument, which has to provide the
  // following methods:
  // - model.f(t,y)
  // - model.df(t,y)
  //
  // replace this function by a function that works for an arbitrary
  // DIRK method
  template <class Model> // one step forward
  double evolve(const double &y, double t, double h, const Model &model) const
  {
    double ret = y; // u_n + h*F in Runge Kutta, set as y becuase we start at u_n
                    // ret is what will be returned
    std::vector<double> k(stages_);
    for (int i = 0; i < stages_; ++i)
    {
      double temp_sum = y; // is the sum u_n + sum_{j=1}^{i-1}a_{ij}K_j
      for (int j = 0; j < i; ++j)
      {
        temp_sum += h*a(i,j)*k[j];
      }
      if (a(i,i) != 0) // if this holds true, implement the newton method.
      {
        int iter = 0;
        double error = -k[i] + model.f(t + h*c_[i], temp_sum + h*a(i,i)*k[i]);
        while (iter < 1e6 && std::abs(error) > 1e-6) // 1e6 = maximum number of iterations and 1e-6 = tolerance for the computation
        {
          double denominator = -1 + model.df(t + h*c_[i], temp_sum + h*a(i,i)*k[i])*h*a(i,i);
          k[i] = k[i] - error/denominator;
          error = -k[i] + model.f(t + h*c_[i], temp_sum + h*a(i,i)*k[i]);
          iter++;
        }
        temp_sum += h*a(i,i)*k[i];
      }
      k[i] = model.f(t+c_[i]*h, temp_sum);
      ret += h*b_[i]*k[i];
    }
    return ret; // return u_{n+1}
  }

protected:

  // accessing A for convenience
  double a(int i, int j) const { return a_[i*stages_+j]; }
  double& a(int i, int j) { return a_[i*stages_+j]; }

  int stages_;

  std::vector<double> a_,b_,c_;
};

// FE: derivied class for the forward Euler method
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

// BE: derived class for the backward Euler method
class BE: public DIRK
{
public:

  BE() : DIRK(1)
  {
	a(0,0) = 1.;
	b_[0] = 1.;
	c_[0] = 1.;
  }

};

class IM: public DIRK
{
public:

  IM() : DIRK(1)
  {
	a(0,0) = 0.5;
	b_[0] = 1.;
	c_[0] = 0.5;
  }

};

class Heun3: public DIRK
{
public:

  Heun3() : DIRK(3)
  {
	a(1,0) = 1./3;
  a(2,1) = 2./3;
	b_[0] = 0.25;
  b_[1] = 0.;
  b_[2] = 0.75;
	c_[0] = 0;
  c_[1] = 1./3;
  c_[2] = 2./3;
  }

};

class DIRK2: public DIRK
{
public:

  DIRK2() : DIRK(3)
  {
  double  delta = 0.5 + sqrt(3)/6;
	a(0,0) = delta;
  a(1,0) = 1 - 2*delta;
  a(1,1) = delta;
	b_[0] = (0.5 - delta)/(1 - 2*delta);
  b_[1] = (0.5 - delta)/(1 - 2*delta);
	c_[0] = delta;
  c_[1] = 1 - delta;
  }

};

#endif
