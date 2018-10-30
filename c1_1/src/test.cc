
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



Sparse m(3,4);
m.addEntry(1,1,0);

std::cout << m.getEntry(1,0) << std::endl;
std::cout << " " << std::endl;

m.checkDiagonal();
//m.checkDimension();


m.printMatrix();

}
