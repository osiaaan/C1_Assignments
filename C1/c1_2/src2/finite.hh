
/*
 file: finite.hh
*/

 #ifndef CLASS_FINITE
 #define CLASS_FINITE

 #include<iostream>
 #include<cmath>
 #include<vector>
 #include<math.h>
 #include<iomanip>
 #include<fstream>
 #include<chrono>
 #include<tgmath.h>
 #include"sparse.hh"
 #include"vector.hh"

 /*
 class: CLASS_FINITE
a class to create a matrix representing the advetion diffusion reaction equation problem
in matrix format along with other fucntionalities
 */

 class Finite
 {
 public:

   //CONSTRUCTORS

   Finite(); // default constructor
   Finite(const Finite& source); //copy constructor;
   Finite(double a, double b, double c, int N); //custom constructor

   //GETTERS
   Sparse getMatrix();
   std::vector<double> getVector();
   std::vector<double> getSolution();

   //FUNCTIONS
   std::vector<double> solve(std::vector<double> u);


private:
   /*
   The Advection Diffusion Reaction equation we are dealing with is of the form
   -au" + bu' + cu = 0.
   We will describe the finite difference method in compact matrix form AU = f.
   */

   double a_, b_ ,c_;//These are the coefficients in the equation
   /*
    The followng is a vector containing the values of the analytical solution to
    the linear advection diffusion equation at the grid points.
   */
   std::vector<double> analyticSolution_;
   int N_;//This will be the mesh size, specifically interior points
   std::vector<double> f_;//The is corresponds the the vector f described in the comment
   Sparse A_;//This corresponds to the matrix A described in the comment

 };

  double analyticSolution(double x, double p);
  double error(std::vector<double> u, std::vector<double> v);
  double test(int N, double beta);

 #endif
