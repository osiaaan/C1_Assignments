#include<iostream>
#include<cmath>
#include<vector>
#include<math.h>
#include<iomanip>
#include<fstream>
#include<chrono>
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
Finite x(0,2,2,1000);
Sparse A = x.getMatrix();
std::vector<double> f = x.getVector();
std::vector<double> u(f.size());

for( double n : u )
{
  n = 0.0;
}

A.GaussSeidel(u,f);
std::cout << " " << std::endl;

}
