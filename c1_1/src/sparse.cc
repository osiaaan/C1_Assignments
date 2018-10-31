#include "sparse.hh"

  double const TOL = 1e-3;

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

//This will return the ijth entry of the Matrix
double Sparse::getEntry(int i, int j)
{
  double a = 0.0;
  std::vector<int> index = (*this).getIndex(i); //ith index#
  std::vector<double> row = (*this).getRow(i); //ith row
  int k_0 = 0;//indicator

  if((*this).checkIndex(i) == true)
  {
      for(int k = 0 ; k < (*this).getWidth(); k++)
      {
        if(j == index[k_0])
        {
          a = row[k_0];
          break;
        }
        else if((k == index[k_0]) && (k_0 < index.size() - 1)){k_0++;} //should probably have used an iterator...
      }
  }

  return a;
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

//This will check if we have ith row all zeroes, returning false if so
bool Sparse::checkIndex(int i) //will check row "i+1" due to indexing starting at zero.
{
  bool indicator = true;

  if(i > (*this).getLength() - 1)
  {
    std::cout << "You are checking for a row which does not exist." <<  std::endl;
    std::cout << "The boolean will be returned as false." << std::endl;
    indicator = false;
  }
  else
  {
    std::vector<int> index = (*this).getIndex(i);
    if(index.size() == 0)
    {
      indicator = false;
      std::cout << "This matrix has a zero row. "<< std::endl;
    }
  }

  return indicator;
}

//This will check if our matrix has a zero row, returning false if so
bool Sparse::checkMatrix()
{
  bool indicator = true;
  for(int i = 0 ; i < (*this).getLength(); ++i)
  {
    indicator = (*this).checkIndex(i);
    if(indicator == false)
    {
      break;
    }
  }
}

//This will check if our matrix has a zero in its diagonal, returning false if so
bool Sparse::checkDiagonal()
{
  bool indicator = true;

  for(int i = 0 ; i < (*this).getLength() ; ++i)
  {
    if((*this).checkDimension() == false)
    {
      std::cout << "boolean returned as false." << std::endl;
      indicator = false;
      break;
    }

    double a = (*this).getEntry(i,i);
    if(a == 0)
    {
      std::cout << "There is a zero in the diagonal of this matrix." << std::endl;
      indicator = false;
      break;
    }

    if(indicator == true)
    {
      std::cout << "This matrix has non-zero diagonal elements." << std::endl;
    }
    return indicator;
  }
  /*
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
  */
}

//This will check to see if the matrix is nxn, returning false if not
bool Sparse::checkDimension()
{
  bool indicator = true;

  int a((*this).getLength()), b((*this).getWidth());

  if(a != b)
  {
    std::cout<< "This is not an nxn matrix!" << std::endl;
    indicator = false;
  }
  else
  {
    std::cout <<"This is an nxn matrix!" <<std::endl;
  }

  return indicator;
}

//This will compute Ax
std::vector<double> Sparse::matrixVectorMult(std::vector<double>x)
{
int n = (*this).getLength();
int m = (*this).getWidth();
std::vector<double> y (n);
double sum;

for(int i = 0 ; i < n ; ++i)
{
  sum = 0.0;
  for(int j = 0 ; j < m ; ++j)
  {
    double a = (*this).getEntry(i,j);
    sum += x[j]*a;
  }
  y[i] = sum;
}

return y;
}

//This will help us calculate the error later on
double infinityNorm(std::vector<double> x) //take a wild guess
{
  double a = abs(x[0]);
  for(int i = 1; i < x.size() ; ++i)
  {
      if(abs(x[i]) > a)
      {
        a = abs(x[i]);
      }
  }
  return a;
}

//This will let us negate vector y from vector x
std::vector<double> minus(std::vector<double> x, std::vector<double> y)
{
  if(x.size() != y.size())
  {
    std::cout << "These vectors are not of the same size. Returned original vector. \n";
    return x;
  }
  else
  {
  for( int i = 0 ; i < x.size() ; ++i )
    {
      x[i] = x[i] - y[i];
    }
    return x;
  }
}

//Here we implement out Gauss-Seidel algorithm.


//std::vector<double> Sparse::GaussSeidel(std::vector<double> x_0, std::vector<double> b)
//{
  /*
  bool diagonal_indicator = (*this).checkDiagonal();
  bool dimensionIndicator = (*this).checkDimension();

  //Here we put a while loop for tolerance
  int n = (*this).getLength();
  int MaxIter = 1000;
  int iter = 0;

return x_0;
 */

/*
 std::vector<double> x_k=  x_0; //our approximation
 std::vector<double> r = minus(b,(*this).matrixVectorMult(x_k)); //The residaul
 int n = (*this).getLength();
 int it = 0;

 for(int l = 0; l < 1000; l++)
 {
   for(unsigned int i = 0; i < n; ++i)
   {
     double sum = 0.0;
     for(unsigned int j = 0; j < n ; j++ )
     {
       if(j != i)
       {
       sum += x_k[j]*( (*this).getEntry(i,j) );
       }
     }
     x_k[i] = (b[i] - sum)/((*this).getEntry(i,i));
   }
   r = minus(b, (*this).matrixVectorMult(x_k));
   it++;

 } //while(infinityNorm(r) > TOL);
 std::cout << "number of iterations: " << it << std::endl;

 return x_k;
}

*/
