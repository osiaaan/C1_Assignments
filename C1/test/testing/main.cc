#include <iostream>
#include "sparse.hh"
#include "vector.hh"

int main()
{
  Vector l(10,2.55);
  std::cout << l[1] << std::endl;

  SparseMatrix A(5);

  return 0;
}
