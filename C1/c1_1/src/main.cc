
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

std::vector<double> delta = {1};

for(int i = 0 ; i < delta.size() ; ++i)
{
  const unsigned int N = 100;//dimension
  const  double delta_ = delta[i];//constant
  double a = 4*(1-delta_);

  Sparse A(N,N);//This will be our tridiagonal matrix

  std::vector<double> w (N);//vector w
  std::vector<double> x (N);//This will be our initial guess
  std::vector<double> b (N);//Ax=b

  //setting up our vector w
  for(int i = 0; i < N; ++i)
  {
    w[i] = (double) (i+1)/(N+1); //casting as double
  }

  //setting up our tridiagonal matrix
  std::vector<double> D (N+1); //plus one to include "D_-1"
  double d = (w[0] - (1.0/2.0));
  D[0] = a*d*d + delta_;//setting up "D_-1 = D_0"

  for(int i = 1 ; i < N + 1 ; i++)//vector gives us the D_is
  {
      double d = (w[i-1] - (1.0/2.0));
      D[i] = a*d*d + delta_;
    }

  for(int i=0 ; i < N ; ++i)//adding the entries to our matrix
  {
      for(int j = 0; j < N; ++j)
      {
        if(j == i-1)
        {
          A.addEntry(-D[i],i,j);
        }
        else if(j == i)
        {
          A.addEntry(D[i+1]+D[i],i,j);
        }
        else if(j == i +1)
        {
          A.addEntry(-D[i+1],i,j);
        }
      }
    }

  //setting up our vector b
  for(int i = 0; i < N -1 ; ++i)
  {
    double c = (w[i] - (1.0/2.0));
    b[i] = (-2)*a*c*w[0]*w[0];
  }
  b[N-1] = (-2)*a*(w[N-1] - (1.0/2.0))*w[0]*w[0] + 1;


  //setting our initial guess to be the zero vector
  for( double n : x )
  {
    n= 0.0;
  }

  /*Here we implement Gauss Seidel, and it returns a vector of the
  residual error for each iteration */
  std::vector<std::vector<double>> result = A.GaussSeidel(x,b);
  std::vector<double> x_sol = result[0];
  std::vector<double> residual = result[1];
  std::cout << "The error between the solution and the approximation is: " << infinityNorm(minus(w,x_sol)) << std::endl;
}


}
