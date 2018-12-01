#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>
#include "models.hh"

double sum(std::vector<double> v)
{
  double s = 0.;
  for(int i = 0 ; i < v.size() ; ++i)
  {
    s +=v [i];
  }
  return s;
}

template <class Model>
double newton (double t, const double &y, double h, const Model &model, std::vector<double> a, std::vector<double> summand, int i, int stages)
{
  double c_0 = 0.;

  double const tol = 1e-6;
	for( int i = 0 ; i < 10000 ; i++ )
	{
		double x = c_0;
		double fx = model.f(t,x) - ( x - h*sum(summand) - y ) / ( h*a[i*stages + i] );
		double fx1 = model.df(t,x) - 1/(h*a[i*stages + i]);
		c_0 = x - (fx/fx1);


		if(std::abs(c_0 - x) < tol)
		{
			break;
		}
	}
  //std::cout << c_0 << std::endl;
	return c_0;
}
