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
  int const N = 100;
//  Finite x(1,1,0,N);

  //Vector u = x.solve();
  //Vector u_sol = x.getSolution();

  //std::cout << "There error between the analytic and approximate solution is: " << infinityNorm(minus(u_sol,u)) << std::endl;
  //data(dat,"solution");

}
