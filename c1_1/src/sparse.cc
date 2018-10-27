#include "sparse.hh"

//constructors

Sparse::Sparse() : N_length_(3), M_width_(3)
{
  // By default we get the 3x3 identity matrix
std::vector<double> v1(1),v2(1),v3(1);
std::vector<int> i1(3),i2(3),i3(3);
v1 = {1};
v2 = {1};
v3 = {1};
i1 = {0};
i2 = {1};
i3 = {2};
sparse_matrix_ = {v1,v2,v3};
indexing_vector_ = {i1,i2,i3};
}

//This is the copy constructor
Sparse::Sparse(const Sparse& source) : sparse_matrix_(source.sparse_matrix_), indexing_vector_(source.indexing_vector_), N_length_(source.N_length_), M_width_(source.M_width_)
{}

//This is the custom constructor
/*
Sparse::Sparse(int N, int M, std::vector<std::vector<double>> sparseRows, std::vector<std::vector<int>> index) : sparse_matrix_(sparseRows), indexing_vector_(index)
{
N_length_ = N;
M_width_ = M;
}
*/
Sparse::Sparse(const int N,const  int M) : N_length_(N), M_width_(M)
{
 sparse_matrix_.resize(N);
 indexing_vector_.resize(N);
}


int Sparse::getLength()
{
  int a = N_length_;
  return a;
}

int Sparse::getWidth()
{
  int a = M_width_;
  return a;
}

std::vector<double> Sparse::getRow(int i)
{
  std::vector<double> v;
  v = sparse_matrix_[i];
  return v;
}

std::vector<int> Sparse::getIndex(int i)
{
  std::vector<int> v;
  v = indexing_vector_[i];
  return v;
}


void Sparse::addEntry(double a, int i, int j)
{
  std::vector<int> u = (*this).getIndex(i);
  std::vector<double> r = (*this).getRow(i);
  sparse_matrix_[i] = r;
  indexing_vector_[i] = u;
}


void Sparse::printMatrix()
{
  for(int i = 0 ; i < (*this).getLength() ; ++i) // choosing row i
  {
    std::vector<double> row_i = (*this).getRow(i); //setting row_i to be the ith vector of sparse_matrix_
    std::vector<int> index_i = (*this).getIndex(i);//setting index_i to be the ith vector of indexing_vector_
    int k = 0;//this will help us keep track of whether to print out a zero or an entry from sparse_matrix_

    if(index_i.size() == 0)
    {
      for(int j = 0 ; j < (*this).getWidth() ; ++j)
      {
        std::cout<< "0 ";
      }
      std::cout << " " << std::endl;
      continue;
    }
    else
    {

    for(int j = 0 ; j < (*this).getWidth() ; ++j)
    {
      if(j == index_i[k])
      {
        std::cout << row_i[k] << " ";
        k += 1;//now looking at the next entry of index_i
      }
      else
      {
        std::cout << "0 ";
      }
    }
    std::cout << " " << std::endl;
    }
  }


}
