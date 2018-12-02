#ifndef SCHEMES_HH
#define SCHEMES_HH

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>
#include "newton.hh"

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
  //
  // replace this function by a function that works for an arbitrary
  // DIRK method
  template <class Model>
  double evolve(const double &y,double t,double h,const Model &model) const
  {
    double ret = y;
    std::vector<double> k(stages_);

    //ignore what is commented out here
    /*
    k[0] = model.f(t+c_[0]*h, y);
    ret += h*b_[0]*k[0];
    return ret;

    double newt = newton(t + c_[0]*h, y, h, model, a_, k, 0, stages_);
    double K_1 = (newt - y)/h;
    k[0] = K_1;
    */

    //Check if the shceme is implicit or not (i.e. FE or not)
    //Note: there is a new protected boolean called implicit
    if(implicit_ == false)
    {
      k[0] = model.f(t+c_[0]*h, y);
      ret += h*b_[0]*k[0];
      return ret;
    }
    else //The shceme is implicit...
    {
      for(int i = 0 ; i < stages_ ; ++i)
      {
        //This will be the value of sum_{j=1}^{i-1}a_{ij}K_j
        double sum;
        //This vector will hold the entries of the sum sum_{j=1}^{i-1}a_{ij}K_j
        std::vector<double> summand (i);
        for(int j = 0 ; j < summand.size() ; ++j)
        {
            summand[j] = a(i,j)*k[j];
        }
        //Using the sum_up function which is defined in newton.hh
        sum = sum_up(summand);

        //We can now use these parameters to solve for K_i via newton
        double newt = newton(t + c_[i]*h, y, h, model, a(i,i), sum);
        double K_i = ( newt - h*sum - y ) / ( h*a(i,i) );
        //Finally set plug this K_i value into our k vector for safekeeping...
        k[i] = K_i;
      }
      for(int i = 0 ; i < stages_ ; ++i)
      {
        // We now return the value of u_{n+1} according to the RK scheme...
        return ret += h*b_[i]*k[i];
      }
    }

  }

protected:

  // accessing A for convenience
  double a(int i, int j) const { return a_[i*stages_+j]; }
  double& a(int i, int j) { return a_[i*stages_+j]; }

  bool implicit_;
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
  implicit_ = false;
  }

};

class BE: public DIRK
{
public:

  BE() : DIRK(1)
  {
	a(0,0) = 1.;
	b_[0] = 1.;
	c_[0] = 1.;
  implicit_= true;
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
  implicit_= true;
  }
};

class Heun3: public DIRK
{
public:

  Heun3() : DIRK(3)
  {
	a(1,0) = 1./3.;
  a(2,1) = 2./3.;
	b_[0] = 0.25;
  b_[1] = 0.;
  b_[2] = 0.75;
  c_[0] = 0,;
  c_[1] = 1./3.;
  c_[2] = 2./3.;
  implicit_= true;
  }
};

class DIRK2: public DIRK
{
public:

  DIRK2() : DIRK(2)
  {
  double delta = 0.5 + (sqrt(3)/6);
	a(0,0) = delta;
  a(1,0) = 1 - 2*delta;
  a(1,1) = delta;
	b_[0] = (0.5 - delta)/(1 - 2*delta);
  b_[0] = (0.5 - delta)/(1 - 2*delta);
  c_[0] = delta;
  c_[1] = 1 - delta;
  implicit_= true;
  }
};

#endif
