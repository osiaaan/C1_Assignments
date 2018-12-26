#include <iostream>
#include <fstream>
#include <sstream>

#include "models.hh"
#include "schemes.hh"

template <class Model>
Vector solve(const Model &model, const DIRK &scheme, double tau)
{
  Vector y=model.y0();

  double t=0.0;
  while ( std::abs(model.T()-t)>tau )
  {
    y =  scheme.evolve(y,t,tau,model);
    t += tau;
    std::cout << t << std::endl;
  }
  std::cout << "Finished time integration at t=" << t << std::endl;

  return y;

}

int main(int argc, char* argv[])
{

  // read parameters as command line arguments
  int N = 16;
  double kappa = 1.0, tau = 0.0005;
  HeatEquation model(N, kappa);
  Heun3 scheme;

  Vector solver(N);
  solver = solve(model, scheme, tau);

  for(unsigned int i = 0; i < N ; ++i)
  {
    std::cout << solver[i] << std::endl;
  }
  // example invocation of solver
  //solve(model,*scheme,tau).toFile("out.dat");

  return 0;

}
