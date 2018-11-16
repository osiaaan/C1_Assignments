#include "sparse.hh"

//This is the fixed tolerance for the ierations in Gauss Seidel
double const TOL = 1e-6;

//CONSTRUCTORS

//This is the default constructor. It gives us the 3x3 identity.
Sparse::Sparse() : N_length_(3), M_width_(3)
{
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
  double a = 0.0; //the default returned value will be zero
  std::vector<int> index = indexing_vector_[i]; //ith index
  std::vector<double> row = sparse_matrix_[i]; //ith row
  std::vector<int>::iterator it; //iterator for index
  std::vector<double>::iterator ro = row.begin(); //iterator for row

  //first we check if the row has any non-zero values
  if((*this).checkIndex(i) == true)
  {
      for(it=index.begin(); it!=index.end() ;++it)
      {
        //check if we've found the jth element of the ith row,...
        if(j == *it)
        {
          //...if so, we return that value
          a = *ro;
          break;
        }
        else{++ro;} //if necessary, we move to the next element of the row vector
      }
  }

  return a;
}

//OPERATOR OVERLOADING

//This will let us compute Ax
std::vector<double> Sparse::operator*(std::vector<double> x)
{
  if(x.size() != (*this).getWidth())
  {
    std::cout << "Dimensions for multiplication are not correct." << std::endl;
    std::cout << "Original vector returned." << std::endl;
    return x;
  }

  int n = N_length_;
  int m = M_width_;
  std::vector<double> y (n);
  double sum;

  for(int i = 0 ; i < n ; ++i)
  {
    sum = 0.0;
    for(int j = 0 ; j < sparse_matrix_[i].size()  ; ++j)
    {
      double a = sparse_matrix_[i][j];
      sum += x[indexing_vector_[i][j]]*a;
    }
    y[i] = sum;
  }
  return y;
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

  for(int i = 0 ; i < N_length_ ; ++i) // choosing row i
  {
    std::vector<double> row_i = sparse_matrix_[i]; //setting row_i to be the ith vector of sparse_matrix_
    std::vector<int> index_i = indexing_vector_[i];//setting index_i to be the ith vector of indexing_vector_
    std::vector<int>::iterator it = index_i.begin(); //iterator for index
    std::vector<double>::iterator ro = row_i.begin();//iterator for the row

    //If the we have a zero row then we print all zeroes
    if(index_i.size() == 0)
    {
      for(int j = 0 ; j < M_width_ ; ++j)
      {
        std::cout << std::setw(3);
        std::cout << "0";
      }
      std::cout << " " << std::endl;
      continue;
    }
    else
    {
    for(int j = 0 ; j < M_width_ ; ++j)
    {
      if(j == *it)//if there's a value, we print it
      {
        std::cout << std::setw(3);
        std::cout << *ro;
        it++;//now looking at the next entry of index_i
        ro++;//and also move to the next value of row
      }
      else//otherwise we print out a zero
      {
        std::cout << std::setw(3);
        std::cout <<  "0";
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
  if((*this).checkDimension() == false)
  {
    std::cout << "boolean returned as false." << std::endl;
    indicator = false;
    return indicator;
  }

  for(int i = 0 ; i < (*this).getLength() ; ++i)
  {
    double a = (*this).getEntry(i,i);

    if(a == 0)
    {
      std::cout << "There is a zero in the diagonal of this matrix." << std::endl;
      indicator = false;
      break;
    }
  }

    if(indicator == true)
    {
      std::cout << "This matrix has non-zero diagonal elements." << std::endl;
    }
    return indicator;
  }

//This will check to see if the matrix isn't nxn, returning false if so
bool Sparse::checkDimension()
{
  bool indicator = true;

  int a((*this).getLength()), b((*this).getWidth());

  if(a != b)
  {
    std::cout<< "This is not an nxn matrix." << std::endl;
    indicator = false;
  }
  else
  {
    std::cout <<"This is an nxn matrix." <<std::endl;
  }

  return indicator;
}


std::vector<double> Sparse::diag()
{
  std::vector<double> d(N_length_);
  for(int i = 0 ; i < N_length_ ; ++i)
  {
    d[i] = 0;
    for(int j = 0 ; j < i + 1; j++)
    {
      if(indexing_vector_[i][j] == i)
      {
        d[i] = sparse_matrix_[i][j];
        break;
      }
    }
  }
  return d;
}

//This will help us calculate the error later on
double infinityNorm(std::vector<double> x) //take a wild guess
{
  double a = fabs(x[0]);
  for(int i = 1; i < x.size() ; ++i)
  {
      if(fabs(x[i]) > a)
      {
        a = fabs(x[i]);
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

//Prints out our vector
template <class T>
void printVector(std::vector<T> v)
{
  for( T n : v )
  {
    std::cout << n << " ";
  }
  std::cout << '\n';
}

/*This function will make data files for the residual error
with respect to different parameters.
These will be used later by gnuplot.
The string will make the data file name distinct.*/
void data(std::vector<std::vector<double>> R, std::string s)
{

  for(int j = 0 ; j < R.size() ; ++j)
  {
    std::vector<double> y = R[j];//y is the jth vector
    std::vector<double> index;

    //This is indexing the iterations
    for( int i = 0 ; i < y.size() ; i++ )
    {
       index.emplace_back(i);
    }

    /*The index is turned into a string so we can reference
    it in the data file.*/
    std::string myStr = std::to_string(j);
    //This will help us distinguish data files from each other.
    std::string res = "data_" + s + myStr + ".dat";
    std::ofstream myFile;
    //Here we open a file with the desired name
    myFile.open(res.c_str());

    //We check if the file has opened correctly
    if( !myFile.good() )
    {
       std::cout << "Failed to open file." << std::endl;
    }

    /*Here we input the data.
    This is the iteration and its respective residual error.*/
    for(int i = 0 ; i < y.size() ; i++)
    {
    myFile << index[i] << '\t' << y[i] << std::endl;
    }
    myFile.close();
  }
}

/*Here we implement out Gauss-Seidel algorithm.
It will return a vector consisting of the residual error
for each iteration in order.
It's size will be the number of iterations*/
std::vector<std::vector<double>> Sparse::GaussSeidel(std::vector<double> x_k, std::vector<double> b)
{
  std::vector<std::vector<double>> result;
  int MaxIter = 500000; //This is our maximum number of iterations
  int it = 0; //This counts the number of iterations

   std::vector<double> r = minus(b,(*this)*x_k); //The residaul
   int n = N_length_;
   std::vector<double> residual;
   std::vector<double> d = (*this).diag();

   //We first check if ther is a zero in the diagonal, and see if dimensions are correct
   if((*this).checkDiagonal() == false)
   {
     std::cout << "The returned vector contains the residual error between b and the initial guess only." << std::endl;
     return result = {x_k, residual};
   }

   //Now we implement our algorithm

   /*While the error is not under tolerance and the Maximum
   number of iterations has not been met... */
   while(infinityNorm(r) > TOL && it < MaxIter)
   {
     for(unsigned int i = 0; i < n; ++i)
     {
       //initialisng our sum to be zero
       double sum = 0.0;
       //for loop over the length of the ith row vector
       for(unsigned int j = 0; j < sparse_matrix_[i].size() ; j++ )
       {
         /*As long as the indexing does not represent a diagonal entry...*/
         if(indexing_vector_[i][j] != i)
         {
           sum += x_k[indexing_vector_[i][j]]*sparse_matrix_[i][j];
         }
       }
       //here we update our approximation
       x_k[i] = (b[i] - sum)/d[i];
     }
     //updating the residual
     r = minus(b, (*this)*x_k);
     //here we add this iterations residual error to our error vector
     residual.emplace_back(infinityNorm(r));
     //documenting that this iteration is finished
     it++;
  }

 //if the maximum number of iterations has been met, we let the user know
 if(it == MaxIter)
 {
    std::cout << "Maximum number of iterations reached. Approximation not within tolerance." << std::endl;
    std::cout << "number of iterations: " << it << std::endl;
    std::cout << "The residual error is: " << infinityNorm(r) << std::endl;
 }
 //if our residual error is within tolerance...
 else
 {
   std::cout << "number of iterations: " << it << std::endl;
   std::cout << "The residual error is: " << infinityNorm(r) << std::endl;
 }

 result = {x_k, residual};
 return result;
}
