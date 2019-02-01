 #include"finite.hh"

 //CONSTRUCTORS

//For the default construct all coefficients are zero.
//We have a 3x3 zero matrix.
 //Finite::Finite() : N_(3), a_(1.0), b_(0.0), c_(0.0)
/*
 {
   double p = b_/a_;
   double q = c_/a_;
   double h = (double) 1/(N_+1);

   double u(1);
   double v(-2);
   double w(1);

   Sparse A(3,3);
   std::vector<double> f(3);
   f[2] = -1;

   std::vector<double> u_sol(N_);
   for( int i = 0; i < N_; ++i )
   {
     u_sol[i] = analyticSolution(i*h,p);
   }
   u_sol.emplace_back(1);

   for(int i = 0; i < N_; ++i)
   {
     A.addEntry(v, i, i);
   }

   for(int i = 0; i < N_-1; ++i)
   {
     A.addEntry(u, i, i+1);
     A.addEntry(w, i+1, i);
     f[i] = 0;
   }

   analyticSolution_ = u_sol;
   f_ = f;
   A_ = A;
 }
 */

 //This is the copy constructor
 Finite::Finite(const Finite& source) : f_(source.f_), A_(source.A_), a_(source.a_), b_(source.b_), c_(source.c_), N_(source.N_)
 {}

 //This is the custom constructor
 Finite::Finite(double a, double b, double c, int N) : N_(N), a_(a), b_(b), c_(c)
 {

   if ( a == 0 )
   {
     std::cout << "The coefficent 'a' should be nonzero." << std::endl;
     std::cout << "The case for a = 1, b = 0, c = 0 has been returned." << std::endl;
     Finite x(1,0,0,N);
     A_ = x.getMatrix();
     f_ = x.getVector();
   }
   else
   {
     double h = (double) 1 / (N_+1);
     double p = b_/a_;
     double q = c_/a_;

     double lo = -(a_/(h*h)) - (b_/(2*h));
     double diag = (2*a_) / (h*h) + c_;
     double up = -(a_ / ( h*h )) + (b_ / ( 2*h ));

     SparseMatrix A(N_);

     Vector f(N,0);
     f[N_ - 1] = -up;

     Vector u_sol(N_+2,0);
     for( int i = 0; i < N_+2; ++i )
     {
       u_sol[i] = analyticSolution(i*h,p);
     }

     for(int i = 0; i < N_; ++i)
     {
       A.addEntry(diag, i, i);
     }

     for(int i = 0; i < N_-1; ++i)
     {
       A.addEntry(lo, i, i+1);
       A.addEntry(up, i+1, i);
     }

     analyticSolution_ = u_sol;
     f_ = f;
     A_ = A;
    }

 }

 //GETTERS
 SparseMatrix Finite::getMatrix()
 {
  return A_;
 }

 Vector Finite::getVector()
 {
   return f_;
 }

 Vector Finite::getSolution()
 {
   return analyticSolution_;
 }

 //FUNCTIONS
 double analyticSolution(double x, double p)
 {
   double e = exp(p);
   double A = 1/(e-1);

   return A*(exp(p*x)-1);
 }

 Vector Finite::solve()
 {
  Vector u(N_,0);
  return (A_.GaussSeidel(f_, u, 1e-6, 100000));
 }
