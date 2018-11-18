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

     double up = -(a_/(h*h)) - (b_/(2*h));
     double diag = (2*a_) / (h*h) + c_;
     double lo = -(a_ / ( h*h )) + (b_ / ( 2*h ));

     Sparse A(N_,N_);

     std::vector<double> f(N_,0);
     f[N_ -1] = -lo;

     std::vector<double> u_sol(N_+2);
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
 Sparse Finite::getMatrix()
 {
   return A_;
 }

 std::vector<double> Finite::getVector()
 {
   return f_;
 }

 std::vector<double> Finite::getSolution()
 {
  return analyticSolution_;
 }

 //FUNCTIONS
 double analyticSolution(double x, double p)
 {
   if(p == 0)
   {
     return x;
   }
   else
   {
     double e = exp(p);
     double A = 1/(e-1);

     return A*(exp(p*x)-1);
   }
 }

 double error(std::vector<double> u, std::vector<double> v)
 {
   return infinityNorm(minus(u,v));
 }

 void test(int N, double beta)
 {
   std::cout << "N: " << N << " Beta/Alpha: " << beta << std::endl;
   int beta_ = floor(beta);
   std::string myStr1 = std::to_string(beta);
   std::string myStr2 = std::to_string(N);
   std::string name = "solution_" + myStr1 + "_" + myStr2;

   Finite x(1,beta,0,N);

   std::vector<double> u_sol = x.getSolution();

   std::vector<double> uGuess(N);
   for( double n : uGuess )
   {
     n = 0.0;
   }

   std::vector<double>  u = x.solve(uGuess);

   std::vector<std::vector<double>> dat = {u_sol, u};

   std::cout << "The error between the analytic and approximate solution is: " << error(u_sol,u) << std::endl;
   data(dat, name);

     std::ofstream myFile;
     myFile.open("error_.dat",std::ios::app);

     //We check if the file has opened correctly
     if( !myFile.good() )
     {
        std::cout << "Failed to open file." << std::endl;
     }

     /*Here we input the data.
     This is the iteration and its respective residual error.*/

     myFile << N << '\t' << beta << '\t' << error(u_sol,u) <<  std::endl;

     myFile.close();

   std::cout << " " << std::endl;
}

 std::vector<double> Finite::solve(std::vector<double> u)
 {
   std::vector<std::vector<double>> V = A_.GaussSeidel(u,f_);
   V[0].insert(V[0].begin(), 0);
   V[0].emplace_back(1);

   return V[0];
 }
