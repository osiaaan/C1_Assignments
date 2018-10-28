
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

Sparse m(5,6);
m.addEntry(4.0,1,0);
m.addEntry(4.0,1,2);
m.addEntry(3.0,1,1);
m.addEntry(2.0,1,0);

m.printMatrix();

}
