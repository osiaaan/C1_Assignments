#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

#include "models.hh"
#include "schemes.hh"

template <class Model>
std::vector<Vector> solve(const Model &model, const DIRK &scheme, double tau)
{
  std::vector<Vector> Y;
  Vector y = model.y0();

  double t=0.0;

  while ( std::abs(model.T()-t)>tau )
  {
    Y.emplace_back(y);
    y =  scheme.evolve(y,t,tau,model);
    t += tau;
  }
  std::cout << "Finished time integration at t=" << t << std::endl;

  return Y;
}

/*
void testing(int sch, double tau, int N, double kappa)
{
  //Ensure the user picks a valid option
  if(( sch < 1 || sch> 3))
  {
    std::cout << "Please choose a Scheme between 1 and 3" << std::endl;
  }
  else
  {
    std::string schemeNum = std::to_string(sch);
    std::string myStr = "table_" + schemeNum + ".txt";

    std::ofstream myFile;
    myFile.open(myStr.c_str());

    const int numWidth = 20;
    myFile << std::left << std::setw(numWidth) << "N" << std::left << std::setw(numWidth)  << "Estimate" <<  std::endl;

    HeatEquation model(N, kappa);

    FE sch1;
    BE sch2;
    Heun3 sch3;
    std::vector<DIRK> Sch = {sch1,sch2,sch3};

    Vector solveIt(N);
    solveIt = solve(model, Sch[sch - 1], tau);
  }
}
*/

int main(int argc, char* argv[])
{
  if (argc<4)
  {
    std::cerr << "Usage: " << std::endl
			  << "  " << argv[0] << " <scheme> <tau_0> <N> <kappa>" << std::endl << std::endl
			  << "where" << std::endl
			  << "  <scheme>   scheme index (int)" << std::endl
			  << "  <tau_0>    initial stepsize (double)" << std::endl
			  << "  <N>        spatial stepsize (int)" << std::endl
        << "  <kappa>    coefficient (double)" << std::endl;
    return 1;
  }
  // read parameters as command line arguments
  const int schemeNumber = atoi( argv[1] );
  double tau = atof( argv[2] );
  const int N = atoi( argv[3] );
  double kappa = atof( argv[4] );

  //testing(schemeNumber, tau, N, kappa);

  HeatEquation model(N, kappa);

  FE sch1;
  BE sch2;
  Heun3 sch3;

  std::vector<DIRK> sch = {sch1,sch2,sch3};

  std::vector<Vector> solver;
  solver = solve(model, sch[schemeNumber - 1], tau);

  std::ofstream file;
  file.open("dat.dat");
  double t = 0.0;
  for(int i = 0; i < solver.size(); ++i)
  {
    for ( int j = 0; j<solver[i].size(); ++j )
    {
  	   file << t << " " << ((double) j)/(solver[i].size()-1) << " " << solver[i][j] << std::endl;
    }
    t += tau;
  }
  file.close();


  //solver.toFile("out.dat");
  // example invocation of solver
  //solve(model,*scheme,tau).toFile("out.dat");

  return 0;

}
