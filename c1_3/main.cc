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
std::vector<std::vector<double>> solve(const Model &model, const DIRK &scheme, double tau)
{
  //The first element of EOC involves using the time "tau_0"
  std::vector<double> TAU, ERROR, MAXERROR_COUNTER;

  double maxError = 0;
  double y=model.y0();
  double t=0;
  int count = 0;

  while ( t <= model.T() )
  {
    y =  scheme.evolve(y,t,tau,model,count);


    t += tau;
    double error = y - model.exact(t);
    TAU.push_back(t);
    ERROR.push_back(fabs(error));


    //std::cout << y << std::endl;
    maxError = std::max(maxError,std::abs(error));

  }
  MAXERROR_COUNTER.push_back(maxError);
  MAXERROR_COUNTER.push_back(count);


  std::vector<std::vector<double>> vec = {TAU, ERROR, MAXERROR_COUNTER};
  return vec;
}

void testing(int j, int k, int level, double tau)
{
  if( (k > 5 || k < 1) || (j<1 || j> 2))
  {
    std::cout << "Please choose Test 1 or 2 and a Scheme between 1 and 5" << std::endl;
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
      std::vector<std::vector<double>> vec = solve(model, sch[k-1], tau);
      double maxError = vec[2][0];
      std::cout << "tau: " << tau << " " << "error: "<< maxError << " " << "counter: " << vec[2][1] << std::endl;

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

void plots(int j)
{
  Test model;
  model.model_ = j;

  FE sch1;
  BE sch2;
  IM sch3;
  Heun3 sch4;
  DIRK2 sch5;

  std::vector<DIRK> sch = {sch1,sch2,sch3,sch4,sch5};


  double tau1 = 0.1;
  double tau2 = 0.01;

  std::vector<double> tau = {tau1, tau2};
  for( int k = 0; k < tau.size(); ++k)
  {
    std::vector<std::vector<double>> vec1 = solve(model, sch[0], tau[k]);
    std::vector<std::vector<double>> vec2 = solve(model, sch[1], tau[k]);
    std::vector<std::vector<double>> vec3 = solve(model, sch[2], tau[k]);
    std::vector<std::vector<double>> vec4 = solve(model, sch[3], tau[k]);
    std::vector<std::vector<double>> vec5 = solve(model, sch[4], tau[k]);
    //std::vector<std::vector<double>> vec1 = solve(model, sch[k-1], tau1);

    //std::vector<std::vector<double>> vec2 = solve(model, sch[k-1], tau2);

    std::string t_step = std::to_string(tau[k]);
    std::string num = std::to_string(j);
    //std::string schemeNum = std::to_string(k);
    std::string myStr = "plot_"+ num + "_" + t_step + ".dat";

    std::ofstream myFile;
    myFile.open(myStr.c_str());

    for(int i = 0 ; i < vec1[0].size() ; ++i)
    {
      myFile << vec1[0][i] << '\t' <<  vec1[1][i] << '\t' <<  vec2[1][i] << '\t' <<  vec3[1][i] << '\t' <<  vec4[1][i] << '\t' <<  vec5[1][i] << std::endl;
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
  plots(modelNumber);
  return 0;
}
