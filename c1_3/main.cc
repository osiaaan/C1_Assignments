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
  //The first element of EOC involves using the time "tau_0"

  double maxError = 0;
  double y=model.y0();
  double t=0;

  while ( t <= model.T() )
  {
    y =  scheme.evolve(y,t,tau,model);
    t += tau;
    double error = y - model.exact(t);
    //std::cout << y << std::endl;
    maxError = std::max(maxError,std::abs(error));

  }
  return maxError;
}

void testing(int j, int k, int level, double tau)
{

  FE sch1;
  BE sch2;
  IM sch3;
  Heun3 sch4;
  DIRK2 sch5;

  std::vector<DIRK> sch = {sch1,sch2,sch3,sch4,sch5};

  if(j == 1)
  {
    Test1 model;

    std::vector<double> ERROR(level);
    std::vector<double> EOC(level-1);
    // solve and display error
    for (int i=0;i<level;++i,tau/=2.)
    {
      double maxError = solve(model, sch[k-1], tau);
      std::cout << tau << " " << maxError << std::endl;

      if(i == 0)
      {
        ERROR[i] = maxError;
      }
      else
      {
        ERROR[i] = maxError;
        EOC[i-1] = log(ERROR[i]/ERROR[i-1])/log(tau/(2*tau));
      }
    }
    for (int i=0;i<level-1;++i)
    {
      std::cout << EOC[i] << std::endl;
    }
  }
  else if(j == 2)
  {
    Test2 model;

    std::vector<double> ERROR(level);
    std::vector<double> EOC(level-1);
    // solve and display error
    for (int i=0;i<level;++i,tau/=2.)
    {
      double maxError = solve(model, sch[k-1], tau);
      std::cout << tau << " " << maxError << std::endl;

      if(i == 0)
      {
        ERROR[i] = maxError;
      }
      else
      {
        ERROR[i] = maxError;
        EOC[i-1] = log(ERROR[i]/ERROR[i-1])/log(tau/(2*tau));
      }
    }
    for (int i=0;i<level-1;++i)
    {
      std::cout << EOC[i] << std::endl;
    }
  }
  else {std::cout << "Please choose Test 1 or 2" << std::endl;}


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
			  << "  <J>        number of EOCs (int)" << std::endl;
    return 1;
  }
  const int modelNumber = atoi( argv[1] );
  const int schemeNumber = atoi( argv[2] );
  double tau = atof( argv[3] );
  const int level = atoi( argv[4] );


  // choose correct model and scheme
/*
  FE sch0;
  BE sch1;
  IM sch2;
  Heun3 sch3;
  DIRK2 sch4;

  std::vector<DIRK> sch = {sch0,sch1,sch2,sch3,sch4};

  Test1 model1;
  Test2 model2;

  std::vector<MODELS> = {model1,model2};
  //Heun3 scheme;

  std::vector<double> ERROR(level);
  std::vector<double> EOC(level-1);
  // solve and display error
  for (int i=0;i<level;++i,tau/=2.)
  {
    double maxError = solve(model, sch[schemeNumber], tau);
    std::cout << tau << " " << maxError << std::endl;

    if(i == 0)
    {
      ERROR[i] = maxError;
    }
    else
    {
      ERROR[i] = maxError;
      EOC[i-1] = log(ERROR[i]/ERROR[i-1])/log(tau/(2*tau));
    }
  }

  for (int i=0;i<level-1;++i)
  {
    std::cout << EOC[i] << std::endl;
  }
*/
  testing(modelNumber, schemeNumber, level, tau);
  return 0;
}
