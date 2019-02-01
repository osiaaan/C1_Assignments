#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>
#include "models.hh"

double sum_up(std::vector<double> v)
{
  double s = 0.;
  for(int i = 0 ; i < v.size() ; ++i)
  {
    s +=v [i];
  }
  return s;
}

template <class Model>
double newton (double t, const double &y, double h, const Model &model, double a_ii, double sum, double c_0)
{
  double const tol = 1e-6;
	for( int i = 0 ; i < 1000000 ; i++ )
	{
		double x = c_0;
		double fx = model.f(t,x) - ( x - h*sum - y ) / ( h*a_ii );
		double fx1 = model.df(t,x) - 1/(h*a_ii);
		c_0 = x - (fx/fx1);


		if(std::abs(c_0 - x) < tol)
		{
      return c_0;
			break;
		}

	}
  //std::cout << c_0 << std::endl;
  //std::cout << "max hit" << std::endl;
	return c_0;
}
