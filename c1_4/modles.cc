#include <cassert>

#include "models.hh"

HeatEquation::HeatEquation(int N, double kappa)
{
  y0_ = y0();

  SparseMatrix A(N);

  double h = 1.0/double(N+1);
  double lo_hi = 1*kappa/(h*h);
  double diag = -2*kappa/(h*h);

  for(unsigned int i = 0; i < N ; ++i)
  {
    A.addEntry(i,i,diag);
  }
  for(unsigned int i = 0 ; i < N-1 ; ++i)
  {
    A.addEntry(i,i+1,lo_hi);
    A.addEntry(i+1,i,lo_hi);
  }

  A_ = A;
  N_=N;
}

Vector HeatEquation::f(double t,const Vector &y) const
{
    Vector f = A_*y;
    f = f*t;
    return f;
}

Vector HeatEquation::y0() const
{
  int N = N_;
  double h = 1.0/double(N+1);

  Vector y(N,0);

  for(unsigned int i = 0; i < N ; i++)
  {
    if((i*h > 0.25) && (i*h <= 0.75))
    {
      y[i] = 1;
    }
  }

  return y;
}

const SparseMatrix HeatEquation::df(double t) const
{
  int N = N_;
  SparseMatrix A(N);

  double h = 1.0/double(N+1);
  double lo_hi = 1*t/(h*h);
  double diag = -3*t/(h*h);

  for(unsigned int i = 0; i < N ; ++i)
  {
    A.addEntry(i,i,diag);
  }
  for(unsigned int i = 0 ; i < N-1 ; ++i)
  {
    A.addEntry(i,i+1,lo_hi);
    A.addEntry(i+1,i,lo_hi);
  }
  return A;
}

double HeatEquation::T() const
{
  return 0.1;
}
