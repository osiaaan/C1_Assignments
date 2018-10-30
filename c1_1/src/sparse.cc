#include "sparse.hh"

//CONSTRUCTORS

//This is the default constructor. It gives us the 3x3 identity.
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

//This is the copy constructor.
Sparse::Sparse(const Sparse& source) : sparse_matrix_(source.sparse_matrix_), indexing_vector_(source.indexing_vector_), N_length_(source.N_length_), M_width_(source.M_width_)
{}

//This is the custom constructor.
//Here we give the dimensions of a zero matrix.
Sparse::Sparse(const int N,const  int M) : N_length_(N), M_width_(M)
{
 sparse_matrix_.resize(N);
 indexing_vector_.resize(N);
}

//GETTERS

//This will get the ijth entry.
double Sparse::getEntry(int i, int j)
{
  double a =  0.0;
  return a;  
}

//This returns the length of our matrix.
int Sparse::getLength()
{
  int a = N_length_;
  return a;
}

//This returns the width of our matrix.
int Sparse::getWidth()
{
  int a = M_width_;
  return a;
}

//This returns a vector representing the ith row of our matrix.
std::vector<double> Sparse::getRow(int i)
{
  std::vector<double> v;
  v = sparse_matrix_[i];
  return v;
}

//This returns a vector representing the indexing of the
//ith row of our matrix.
std::vector<int> Sparse::getIndex(int i)
{
  std::vector<int> v;
  v = indexing_vector_[i];
  return v;
}

//FUNCTIONS

//This allows us to put double a into the ijth position of our Matrix.
void Sparse::addEntry(double a, int i, int j)
{
  std::vector<int> u = (*this).getIndex(i); //Here we're getting the index for row i
  std::vector<double> r = (*this).getRow(i);//This is the ith row

  //establishing iterators for index and row resp.
  std::vector<int>::iterator it_u;
  std::vector<double>::iterator it_r;
  //initialisng the iterators to begin at their respective vectors
  it_u = u.begin();
  it_r = r.begin();

  if(u.size() == 0)
  {
    u.insert(it_u,j);
    r.insert(it_r,a);
  }
  else
  {
    for(int l = 0; l < (*this).getWidth(); ++l)
    {
      if(l == j)//if l is the jth column...
      {
        if(l == *it_u)//first check if the index is also equal...
        {
          u.erase(it_u);
          r.erase(it_r);//if so we erase the value so it can be replaced with 'a'
        }

        u.insert(it_u,j);//Here we insert our index
        r.insert(it_r,a);//Here we insert the new value 'a'
      }
      else if(l == *it_u)//if l is equal to the index...
      {
        it_u = it_u + 1;//we move one to the next index
        it_r = it_r +1;//we move to the next nonzero value in the row
      }
    }
  }

  sparse_matrix_[i] = r;//we replace the ith row with the new one
  indexing_vector_[i] = u;//we replace the ith index with the new one
}

//This will print our sparse matrix in a recognisable form.
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

bool Sparse::checkDiagonal()
{
  int n = (*this).getLength();
  int k_0 = 0;
  bool indicator = true;

  for(int i = 0; i < n; ++i)
  {
    std::vector<int> index_i = (*this).getIndex(i);
    std::vector<double> row_i =(*this).getRow(i);
    if(index_i.size() == 0)
    {
      indicator = false;
      break;
    }

    for(int l = 0; l < i ; ++i) //Finding the index corresponding to the a_iith element.
    {
      if(l == index_i[k_0] )
      {
        ++k_0;
      }
    }
    if(i != index_i[k_0])
    {
      indicator = false;
      break;
    }
  }

  if(indicator == false)
  {
    std::cout<< "There was a zero diagonal element" << std::endl;
  }
  return indicator;
}

//Here we implement out Gauss-Seidel algorithm.
/*
std::vector<double> Sparse::GaussSeidel(std::vector<double> x_0, std::vector<double> b)
{

  //First we need to check that the dimensions are correct! TO_DO
  //Also need to do a check so that the diagonal entries are non zero

  //Here we put a while loop for tolerance
int n((*this.getLength())); //Assuming now that A is nxn, x_0 nx1

for(int i = 0; i < n; ++i)
{
  std::vector<int> index_i = (*this).getIndex(i);
  std::vector<double> row_i =(*this).getRow(i);

  double sum = 0.0; //This will be our sum, also ressets to zero
  int k = 0; //This indicator will help us find the nozero elements of the row.
  int k_0 = 0; // This indicator will help us find a_ii index.

  for(int j = 0; j < n; ++j)
  {
  if(j != i)
  {
    if(j == index[k]) //if the ijth element nonzero
    {
    sum += x_0[j]*row_i[k];// past sumand + new entry
    ++k;
    }
  }
  }

  for(int l = 0; l < i+1 ; ++i) //Finding the index corresponding to the a_iith element.
  {
    if(l == index[k_0] )
    {
    ++k_0;
    }
  }
  double oneOver_a_ii = 1/row[k_0];
  x_0[i] = oneOver_a_ii*(b[i] - sum);
}
//Here we end the while loop for tolerance
}
*/
