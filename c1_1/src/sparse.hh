/*
 file: sparse.hh
*/

 #ifndef CLASS_SPARSE
 #define CLASS_SPARSE

 #include<iostream>
 #include<cmath>
 #include<vector>
 #include<math.h>
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
   std::vector<double> operator*(std::vector<double> source);

   //FUNCTIONS

   void printMatrix();
   void addEntry(double a, int i, int j);
   bool checkIndex(int i); //check if we have ith row all zeroes, returns false if so
   bool checkMatrix();//check if it has a zero row, returns false if so
   bool checkDiagonal();//check if matrix has a zero in diagonal, returns false if so
   bool checkDimension();//checks if the matrix is nxn, returns false if not
   std::vector<double> GaussSeidel(std::vector<double> x_0, std::vector<double> b);

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



 #endif
