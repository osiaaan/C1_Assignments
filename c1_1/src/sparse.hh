/*
 file: sparse.hh
*/

 #ifndef CLASS_SPARSE
 #define CLASS_SPARSE

 #include<iostream>
 #include<cmath>
 #include<vector>

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

   //FUNCTIONS

   void printMatrix();
   void addEntry(double a, int i, int j);
   bool checkDiagonal();
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
