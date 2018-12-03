#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include<fstream>
#include <vector>
#include <iomanip>

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
  if( k > 5 || k < 1)
  {
    std::cout << "Please pick a scheme between 1 and 5" << std::endl;
  }
  else
  {

    std::string testNum = std::to_string(j);
    std::string schemeNum = std::to_string(k);
    std::string myStr = "table_" + testNum + "_" + schemeNum;
    std::ofstream myFile;
    myFile.open(myStr.c_str());
    const int numWidth = 20;
    myFile << std::left << std::setw(numWidth) << "Tau" << std::left << std::setw(numWidth)  << "Error" << std::left << std::setw(numWidth)  << "EOC" << std::endl;

    FE sch1;
    BE sch2;
    IM sch3;
    Heun3 sch4;
    DIRK2 sch5;

    std::vector<DIRK> sch = {sch1,sch2,sch3,sch4,sch5};

    Test model;
    model.model_ = j;

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
        myFile << std::left << std::setw(numWidth) << tau << std::left << std::setw(numWidth) << maxError << std::left << std::setw(numWidth)  << " " << std::endl;
      }
      else
      {
        ERROR[i] = maxError;
        EOC[i-1] = log(ERROR[i]/ERROR[i-1])/log(tau/(2*tau));
        myFile << std::left << std::setw(numWidth) << tau << std::left << std::setw(numWidth)  << maxError << std::left << std::setw(numWidth)  << EOC[i-1] << std::endl;
      }
    }
    myFile.close();
  }
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

  testing(modelNumber, schemeNumber, level, tau);
  return 0;
}
