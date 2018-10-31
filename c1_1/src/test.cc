
#include<iostream>
#include<string>
#include<sstream>
#include<math.h>
#include<cmath>
#include<cstdlib>
#include<ctime>
#include<iomanip>
#include<fstream>
#include<vector>
#include"sparse.hh"



void printVector(std::vector<double> v)
{
  for( double n : v )
  {
    std::cout << n << " ";
  }
  std::cout << '\n';
}


 int main()
{

Sparse m(3,3);

std::vector<double> b = {0,0,1};
std::vector<double> x = {0,0,0};

m.addEntry(2,0,0);
m.addEntry(-1,0,1);
m.addEntry(-1,1,0);
m.addEntry(2,1,1);
m.addEntry(-1,1,2);
m.addEntry(-1,2,1);
m.addEntry(2,2,2);

std::vector<double> y = m.GaussSeidel(x,b);

}
