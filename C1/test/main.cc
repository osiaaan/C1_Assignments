#include <cmath>

#include "output.hh"

int main() {
  const int N = 128;

  // allocate some data
  std::vector<double> x(N);
  std::vector<double> f(N);
  std::vector<double> f_prime(N);


  // compute some results
  for (int i=0; i<N; ++i ) {
	x[i] = 2.0*M_PI/N*i;
	f[i] = std::sin(x[i]);
	f_prime[i] = std::cos(x[i]);
  }

  // output some results for post-processing
  outputVector(f, "func.dat");
  outputVector(f_prime, "deriv.dat");

  return 0;
}

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
  }
  std::cout << "Finished time integration at t=" << t << std::endl;

  return y;

}

int main(int argc, char* argv[])
{

  // read parameters as command line arguments
  int N = 16;
  double kappa = 1.0, tau = 0.0001;
  HeatEquation model(N, kappa);
  FE scheme;

  Vector solver(N);
  solver = solve(model, scheme, tau);

  for(unsigned int i = 0; i < N ; ++i)
  {
    std::cout << solver[i] << std::endl;
  }
  solver[0].toFile("ya bish");
  // example invocation of solver
  //solve(model,*scheme,tau).toFile("out.dat");

  return 0;

}
