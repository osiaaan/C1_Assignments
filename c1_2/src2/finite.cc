#include"finite.hh"

 //CONSTRUCTORS

//This is the default constructor. It is identical to the custom constructor but
//specific values of N,a,b,c.
 Finite::Finite() : N_(3), a_(1.0), b_(0.0), c_(0.0)
 {
     //Here we define the necessary variables
     double h = (double) 1/(N_+1);
     double p = b_/a_;

     //These will be the tridiagonal coefficients
     double lo = -(a_/(h*h)) - (b_/(2*h));
     double diag = (2*a_) / (h*h) + c_;
     double up = -(a_ / ( h*h )) + (b_ / ( 2*h ));

     Sparse A(N_,N_);

     std::vector<double> f(N_,0);
     f[N_ -1] = -up;

     //Here we set up the matrix A representing the operator L_h
     for(int i = 0; i < N_; ++i)
     {
       A.addEntry(diag, i, i);
     }

     for(int i = 0; i < N_-1; ++i)
     {
       A.addEntry(up, i, i+1);
       A.addEntry(lo, i+1, i);
     }

     //This vector will contain the analytic solution's value at the gridpoints
     std::vector<double> u_sol(N_+2);
     for( int i = 0; i < N_+2; ++i )
     {
       u_sol[i] = analyticSolution(i*h,p);
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
   //The Peclet number does not make sense for a = 0
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
     //Here we define the necessary variables
     double h = (double) 1/(N_+1);
     double p = b_/a_;

     //These will be the tridiagonal coefficients
     double lo = -(a_/(h*h)) - (b_/(2*h));
     double diag = (2*a_) / (h*h) + c_;
     double up = -(a_ / ( h*h )) + (b_ / ( 2*h ));

     Sparse A(N_,N_);

     std::vector<double> f(N_,0);
     f[N_ -1] = -up;

     //Here we set up the matrix A representing the operator L_h
     for(int i = 0; i < N_; ++i)
     {
       A.addEntry(diag, i, i);
     }

     for(int i = 0; i < N_-1; ++i)
     {
       A.addEntry(up, i, i+1);
       A.addEntry(lo, i+1, i);
     }

     //This vector will contain the analytic solution's value at the gridpoints
     std::vector<double> u_sol(N_+2);
     for( int i = 0; i < N_+2; ++i )
     {
       u_sol[i] = analyticSolution(i*h,p);
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

 //returns the value of the analytic solution at point x
 double analyticSolution(double x, double p)
 {
   //In case p = 0, the coefficient A below is udnefined
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

//This will be the L_infinity norm of the difference of the vectors
 double error(std::vector<double> u, std::vector<double> v)
 {
   return infinityNorm(minus(u,v));
 }

/*The following function is here to automate the creation of the data files.
Moreover, it returns the error between the analytic and approximate solution.

The inputed arguments are the number of interior gridpoints and the value for
beta (alpha is always assumed to be 1 for simplicity). 
The function will create data files for the approximated and analytic solution,
so that they may be plotted via gnuplot.
*/
 double test(int N, double beta)
 {
   //Noting which case we're considering
   std::cout << "N: " << N << " Beta/Alpha: " << beta << std::endl;
   //converting to strings for naming
   std::string myStr1 = std::to_string(beta);
   std::string myStr2 = std::to_string(N);
   std::string name = "solution_" + myStr1 + "_" + myStr2;

   //Initialise the matrix, function etc.
   Finite x(1,beta,0,N);

   //setting up our inital guess, which is the vector of zeroes
   std::vector<double> uGuess(N);
   for( double n : uGuess )
   {
     n = 0.0;
   }

   //The following vecors represent the analytic and approximated
   //solution respectively.
   std::vector<double> u_sol = x.getSolution();
   std::vector<double>  u = x.solve(uGuess);

   //setting up the data files for plotting
   std::vector<std::vector<double>> dat = {u_sol, u};
   //this function was defined in the sparse class and will create
   //the desired dat files for gnuplot.
   data(dat, name);

   //Here we output the error to the terminal
   std::cout << "The error between the analytic and approximate solution is: " << error(u_sol,u) << std::endl;

   std::cout << " " << std::endl;
   return error(u_sol,u);
 }

//This function uses the gauss seidel member function created in the Sparse
//class to solve our system
 std::vector<double> Finite::solve(std::vector<double> u)
 {
   //V[0] is the approximate solution and V[1] is the vector of residual errors
   //Clearly we won't need V[1] here
   std::vector<std::vector<double>> V = A_.GaussSeidel(u,f_);
   //The following is to account for the boundary conditions
   V[0].insert(V[0].begin(), 0);
   V[0].emplace_back(1);

   return V[0];
 }
