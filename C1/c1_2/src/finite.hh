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
 class: CLASS_SPARSE

 a class to create matrices which help described a finite difference method
 for a case of the dvection Diffusion Reaction equation in compact Matrix
 form. This problem will then be solved using GaussSeidel
 */

 class Finite
 {
 public:

   //CONSTRUCTORS

   //Finite(); // default constructor
   Finite(const Finite& source); //copy constructor;
   Finite(double a, double b, double c, int N); //custom constructor

   //GETTERS
   SparseMatrix getMatrix();
   Vector getVector();
   Vector getSolution();

   //OPERATOR OVERLOADING

   //FUNCTIONS
   Vector solve();


private:
   /*
   The Advection Diffusion Reaction equation we are dealing with is of the form
   -au" + bu' + cu = 0.
   We will describe the finite difference method in compact matrix form AU = f.
   */

   double a_, b_ ,c_;//These are the coefficients in the equation
   /*
    The followng is a vector containing the values of the analytical solution to
    the linear advection diffusion equation at the interior grid points.
   */
   Vector analyticSolution_;
   int N_;//This will be the mesh size, specifically interior points
   Vector f_;//The is corresponds the the vector f described in the comment
   SparseMatrix A_;//This corresponds to the matrix A described in the comment

 };

  double analyticSolution(double x, double p);


 #endif
