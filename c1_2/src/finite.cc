 #include"finite.hh"

 //CONSTRUCTORS

//For the default construct all coefficients are zero.
//We have a 3x3 zero matrix.
 Finite::Finite() : N_(3), a_(1.0), b_(0.0), c_(0.0)
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
     double h = (double) 1/(N_+1);
     double p = b_/a_;
     double q = c_/a_;

     double u(1 - (0.5*h*p));
     double v((q*h*h) - 2);
     double w(1 + (0.5*h*p));

     Sparse A(N_,N_);

     std::vector<double> f(N_);
     f[N_ -1] = (0.5*h*p) - 1;

     std::vector<double> u_sol(N_+2);
     for( int i = 0; i < N_+2; ++i )
     {
       u_sol[i] = analyticSolution(i*h,p);
     }

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

 }

 //GETTERS
 Sparse Finite::getMatrix()
 {
   Sparse A = A_;
   return A;
 }

 std::vector<double> Finite::getVector()
 {
   std::vector<double> f = f_;
   return f;
 }

 std::vector<double> Finite::getSolution()
 {
   std::vector<double> s = analyticSolution_;
   return s;
 }

 //FUNCTIONS
 double analyticSolution(double x, double p)
 {
   double e = exp(p);
   double A = 1/(e-1);

   return A*(exp(p*x)-1);
 }

 std::vector<double> Finite::solve(std::vector<double> u, std::vector<double> f)
 {
   std::vector<std::vector<double>> V = A_.GaussSeidel(u,f);
   V[0].insert(V[0].begin(), 0);
   V[0].emplace_back(1);

   return V[0];
 }
