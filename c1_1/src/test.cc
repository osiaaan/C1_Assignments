
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

template <class T>
void printVector(std::vector<T> v)
{
  for( int n : v )
  {
    std::cout << n << " ";
  }
  std::cout << '\n';
}

 int main()
{

Sparse m(3,3);
m.addEntry(4.0,2,2);

m.printMatrix();

}
