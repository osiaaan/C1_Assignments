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
   Sparse(int N, int M);

   //GETTERS

   std::vector<double> getRow(int i);
   std::vector<int> getIndex(int i);
   int getLength();
   int getWidth();
   double getEntry(int i, int j);

   //OPERATOR OVERLOADING

   std::vector<double> operator*(std::vector<double> source); //computes Ax

   //FUNCTIONS

   void printMatrix();
   void addEntry(double a, int i, int j);
   bool checkIndex(int i); //check if we have ith row all zeroes, returns false if so
   bool checkMatrix();//check if it has a zero row, returns false if so
   bool checkDiagonal();//check if matrix has a zero in diagonal, returns false if so
   bool checkDimension();//checks if the matrix is nxn, returns false if not
   std::vector<double> diag();
   std::vector<std::vector<double>> GaussSeidel(std::vector<double> x_k, std::vector<double> b);
  // Sparse tridiagonal(double delta);

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
