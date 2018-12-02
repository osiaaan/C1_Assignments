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
    // putting in forward Euler here temporarily
    /*
    k[0] = model.f(t+c_[0]*h, y);
    ret += h*b_[0]*k[0];
    return ret;
    */
    //double newt = newton(t + c_[0]*h, y, h, model, a_, k, 0, stages_);
    //double K_1 = (newt - y)/h;
    //k[0] = K_1;

    //Check if the shceme is implicit or not (i.e. FE or not)
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
        double sum;
        std::vector<double> summand (i);
        for(int j = 0 ; j < summand.size() ; ++j)
        {
            //This is sum_{j=1}^{i-1}a_{ij}K_j
            summand[j] = a(i,j)*k[j];
        }
        sum = sum_up(summand);

        double newt = newton(t + c_[i]*h, y, h, model, a(i,i), sum);
        double K_i = ( newt - h*sum - y ) / ( h*a(i,i) );
        k[i] = K_i;
      }
      //k[0] = model.f(t + c_[0]*h, y + h*Kay);
      for(int i = 0 ; i < stages_ ; ++i)
      {
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

#endif
