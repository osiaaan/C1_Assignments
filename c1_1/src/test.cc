
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

/*
std::vector<double> v = {1,1,3};
std::vector<double> u = {10,9,3};

Sparse m;
m.addEntry(1,1,0);

std::cout << m.getEntry(1,0) << std::endl;
std::cout << " " << std::endl;

m.checkDiagonal();
//m.checkDimension();



printVector(u);
m.printMatrix();
*/


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
//std::cout << m.getEntry(2,2) << std::endl;
//m.addEntry(2,1,2);
//m.addEntry(1,2,0);
//double a = infinityNorm(b);
//std::cout << a << std::endl;
//if (m.checkIndex(0) == true ){std::cout << "yes" <<std::endl;}
//std::vector<int> y = m.getIndex(0);
//printVector(y);

//Prints Entries of Matrix in order
/*
for(int k = 0; k < m.getLength(); ++k)
{
for(int i = 0; i < m.getWidth(); ++i)
{
std::cout << m.getEntry(k,i) << std::endl;
}
}
*/

std::vector<double> y = m.GaussSeidel(x,b);
m.printMatrix();

std::cout << " " << std::endl;

printVector(y);


}
