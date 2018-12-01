#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>
#include "models.hh"

template <class Model>
double newton (double t, const double &y, double h, const Model &model, double c_0, int maxIter)
{
  double const tol = 1e-6;
	for( int i = 0 ; i < maxIter ; i++ )
	{
		double x = c_0;
		double fx = model.f(t,x) - ((x - y)/h);
		double fx1 = model.df(t,x) - (1/h);
		c_0 = x - (fx/fx1);


		if(std::abs(c_0 - x) < tol)
		{
			break;
		}
	}
  //std::cout << c_0 << std::endl;
	return c_0;
}
/*
for (int i = 0; i < maxIter; i++) {

            double x = c_0;

            double f = model.f(t + c_[0] * h, x) - (x - y) / h;

            double df = model.df(t + c_[0] * h, x) - (1 / h);

            c_0 = x - (f / df);

            if (fabs(c_0 - x) < tol) {
                break;
            }
        }
*/
