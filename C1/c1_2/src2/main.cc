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
  /*The following vectors allow the testing to be automated for
  different mesh size and values of beta (i.e. different Peclet numbers).
  Note that a number of gridpoints 99 gives mesh size h = 1/(99+1) = 0.01 etc.
  */
  std::vector<int> number = {99,499,999};
  std::vector<double> beta = {0, 2 , 20, 200};
  //This vector will collect the required data to be able to plot the error
  //of the approximated solution as mesh size h varies
  std::vector<double> Error (number.size());

  //Here we are just iterating over the vectors and creating dat files for the analytic and approximate solutions
  //to be read by gnuplot
  for(int i = 0; i < beta.size() ; ++i)
  {
    for(int j = 0 ; j < number.size() ; ++j)
    {
      Error[j] = test(number[j],beta[i]);
    }

    //Now we create the dat files for the error vs mesh size
    //to be read by gnuplot

    //setting up a string to distinguish the data files
    int beta_ = floor(beta[i]);
    std::string myStr = std::to_string(beta_);
    std::string myStr2 = "error" + myStr+ ".dat";

    std::ofstream myFile;
    myFile.open(myStr2.c_str());

    if( !myFile.good() )
    {
       std::cout << "Failed to open file." << std::endl;
    }

    //Here we just input the data
    for(int k = 0 ; k < number.size() ; k++)
    {
      double h = (double) 1 /(1 + number[k]);
      myFile << h << '\t' << Error[k] << std::endl;
    }

    myFile.close();
  }

}
