#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

#include "models.hh"
#include "schemes.hh"

template <class Model>
std::vector<Vector> solve(const Model &model, const DIRK &scheme, double tau)
{
  //This vector will contain the approximate solution U at each time-step
  std::vector<Vector> Y;
  //initial condition
  Vector y = model.y0();

  double t=0.0;

  while ( std::abs(model.T()-t)>tau )
  {
    Y.emplace_back(y);
    y =  scheme.evolve(y,t,tau,model);
    t += tau;
  }
  Y.emplace_back(y);

  std::cout << "Finished time integration at t=" << t << std::endl;

  return Y;
}

//This function gives data files which show that the shchemes converge for
//the parameters tau stated in section 2.1.
//With time permiting I would have made tables
void stability(std::vector<Vector> solver, double tau, int schemeNumber)
{
  //These are to distinguish the files
  std::string taau = std::to_string(tau);
  std::string num = std::to_string(schemeNumber);

  std::ofstream file;
  file.open("stabilityData_" + num + "_" + taau +".dat");
  double t = 0.0;

  //We look at whether or not it's blown up at the last time step
  int i = solver.size() -1;

  //This is to account for boundary conditions
  solver[i].push_back(0.0);
  solver[i].emplace_back(0.0);

  for ( int j = 0; j<solver[i].size(); ++j )
    {
       file << ((double) j)/(solver[i].size()-1) << "\t" << solver[i][j] << std::endl;
    }
  file.close();
}


//The following file is to create a 3D plot of the approximated solution
//It's implementation becomes clearer in the scope of main()
void ploting(std::vector<Vector> solver, double tau, int schemeNumber)
{
  std::ofstream file;

  //This is to distinguish the files
  std::string taau = std::to_string(tau);
  std::string num = std::to_string(schemeNumber);
  file.open("plot_" + num + "_" + taau + ".dat");
  double t = 0.0;

  for(int i = 0; i < solver.size(); ++i)
  {
    //These are to account for the initial conditions
    solver[i].push_back(0.0);
    solver[i].emplace_back(0.0);

    for ( int j = 0; j<solver[i].size(); ++j )
    {
       file << t << "\t" << ((double) j)/(solver[i].size()-1) << "\t" << solver[i][j] << std::endl;
    }
    t += tau;
  }
  file.close();
}


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

  //plotting the approximate solution and creating files to "document stability"
  stability(solver, tau, schemeNumber);
  ploting(solver, tau, schemeNumber);
  std::cout << solver[solver.size()-1] << std::endl;


  return 0;

}
