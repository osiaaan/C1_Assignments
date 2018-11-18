#include<iostream>
#include<cmath>
#include<vector>
#include<math.h>
#include<iomanip>
#include<fstream>
#include<chrono>
#include"sparse.hh"
#include"finite.hh"



int main()
{
  std::vector<int> number = {100,1000};
  std::vector<double> beta = {0,2,20,100,200};

  for(int i = 0; i < number.size() ; ++i)
  {
    for(int j = 0 ; j < beta.size() ; ++j)
    {
      test(number[i],beta[j]);
    }
  }

}
