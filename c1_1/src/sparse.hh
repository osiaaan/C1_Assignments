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
   Sparse(); // default constructor
   Sparse(const Sparse& source); //copy constructor;
   Sparse(int N, int M);
   //Sparse(int N, int M, std::vector<std::vector<double>> sparseRows, std::vector<std::vector<int>> index);


   //getters and setters
   //addEntry(double a);

   std::vector<double> getRow(int i);
   std::vector<int> getIndex(int i);
   int getLength();
   int getWidth();

   void printMatrix();
   void addEntry(double a, int i, int j);

private:
   /*
   Here we are setting the two 'vectors of vectors' mentioned in the
   assingment.
   */
   int N_length_;
   int M_width_;
  std::vector<std::vector<double>> sparse_matrix_;
  std::vector<std::vector<int>>  indexing_vector_;

 };

 #endif
