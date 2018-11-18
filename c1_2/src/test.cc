#include<iostream>
#include<cmath>
#include<vector>
#include<math.h>
#include<iomanip>
#include<fstream>
#include<chrono>
#include"sparse.hh"
#include"finite.hh"


template <class T>
void printVector_(std::vector<T> v)
{
  for( T n : v )
  {
    std::cout << n << " ";
  }
  std::cout << '\n';
}

int main()
{
  int const N = 10;
  Finite x(1,1,0,N);

  std::vector<double> u;
  std::vector<double> u_sol = x.getSolution();
  std::vector<double> uGuess(N);

  (x.getMatrix()).printMatrix();

  for( double n : uGuess )
  {
    n = 0.0;
  }

  u = x.solve(uGuess, x.getVector());

  std::vector<std::vector<double>> dat = {u_sol, u};
  std::cout << u.size() << " " << u_sol.size() << std::endl;

  std::cout << "There error between the analytic and approximate solution is: " << infinityNorm(minus(u_sol,u)) << std::endl;
  data(dat,"solution");

}
