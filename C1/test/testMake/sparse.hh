/*
 file: sparse.hh
*/

 #ifndef CLASS_SPARSE
 #define CLASS_SPARSE

 #include<iostream>
 #include<cmath>
 #include<vector>
 #include<math.h>
 #include<iomanip>
 #include<fstream>
 #include<chrono>
 
 /*
 class: CLASS_SPARSE

 a class to hold sparse matricies
 */

 class Sparse
 {
 public:

   //CONSTRUCTORS

   Sparse(); // default constructor
   Sparse(const Sparse& source); //copy constructor;
   Sparse(int N, int M); //custom constructor, gives zero matrix of desired dimesnions

   //GETTERS

   std::vector<double> getRow(int i);
   std::vector<int> getIndex(int i);
   int getLength();
   int getWidth();
   double getEntry(int i, int j);

   //OPERATOR OVERLOADING

   std::vector<double> operator*(std::vector<double> source); //computes Ax

   //FUNCTIONS

   void printMatrix();//prints the matrix in a recognisable form
   void addEntry(double a, int i, int j);//functionality is clear. I suppose this is a setter.
   bool checkIndex(int i); //check if we have ith row all zeroes, returns false if so
   bool checkMatrix();//check if it has a zero row, returns false if so
   bool checkDiagonal();//check if matrix has a zero in diagonal, returns false if so
   bool checkDimension();//checks if the matrix isn't nxn, returns false if so
   std::vector<double> diag();//returns a vector containg the diagonal elements of the matrix
   /*The following implements the Gauss Seidel matrix for guess x_k and vector b.
   It returns a vector of 2 vectors.
   The first is the approximate solution.
   The second is a vector containing the residual error at each iteration. */
   std::vector<std::vector<double>> GaussSeidel(std::vector<double> x_k, std::vector<double> b); //

private:
   /*
   Here we are setting the two 'vectors of vectors' mentioned in the
   assingment, along with the dimensions of the matrix.
   */
   int N_length_;
   int M_width_;
   std::vector<std::vector<double>> sparse_matrix_;
   std::vector<std::vector<int>>  indexing_vector_;
 };

double infinityNorm(std::vector<double> x); //computes the L_infinity norm of the vector

std::vector<double> minus(std::vector<double> x, std::vector<double> y);//negates vector y from x

template <class T>
void printVector(std::vector<T> v);//prints the vector v

void data(std::vector<std::vector<double>> R, std::string s);//creates datafiles residual_"s""i"

 #endif
