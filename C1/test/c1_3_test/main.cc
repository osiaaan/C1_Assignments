#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>

#include "schemes.hh"
#include "models.hh"

template <class Model>
double solve(const Model &model, const DIRK &scheme, double tau)
{
  double maxError = 0;
  double y=model.y0();
  double t=0;
  while ( t<=model.T() )
  {
    y =  scheme.evolve(y,t,tau,model);
    t += tau;
    double error = y - model.exact(t);
    maxError = std::max(maxError,std::abs(error));
  }
  return maxError;
}

int main ( int argc, char **argv )
{
  if (argc<5) 
  {
    std::cerr << "Usage: " << std::endl
			  << "  " << argv[0] << " <model> <scheme> <tau_0> <J>" << std::endl << std::endl
			  << "where" << std::endl
			  << "  <model>    model index (int)" << std::endl
			  << "  <scheme>   scheme index (int)" << std::endl
			  << "  <tau_0>    initial stepsize (double)" << std::endl
			  << "  <J>        number of eocs (int)" << std::endl;
    return 1;
  }
  const int modelNumber = atoi( argv[1] );
  const int schemeNumber = atoi( argv[2] );
  double tau = atof( argv[3] );
  const int level = atoi( argv[4] );


  // choose correct model and scheme
  Test model;
  FE scheme;
  
  // solve and display error
  for (int i=0;i<level;++i,tau/=2.)
  {
    double maxError = solve(model,scheme,tau);
    std::cout << tau << " " << maxError << std::endl;
  }
  return 0;
}
