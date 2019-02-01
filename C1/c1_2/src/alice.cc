#include "FiniteDifferenceSolver.h"


using namespace std;


FiniteDifferenceSolver::FiniteDifferenceSolver() {
}

FiniteDifferenceSolver::FiniteDifferenceSolver(const FiniteDifferenceSolver &source)
        : gridPoints_N_(source.gridPoints_N_), A_(source.A_), alpha_(source.alpha_), beta_(source.beta_),
          gamma_(source.gamma_), L(source.L) {
}

FiniteDifferenceSolver::FiniteDifferenceSolver(int N, double a, double b, double c, double L)
        : gridPoints_N_(N), alpha_(a), beta_(b), gamma_(c), L(L) {

    double h = double(1) / double(N + 1);

    double diag_entry = (2 * alpha_) / (h * h) + gamma_;
    double upper_entry = (-1 * alpha_) / (h * h) + (beta_) / (2 * h);
    double lower_entry = ((-1 * alpha_) / (h * h)) - ((beta_) / (2 * h));

    //cout << lower_entry << endl;

    SparseMatrix A(N);

    //construct the matrix A_

    //add diagonal entries

    for (unsigned int i = 0; i < N; i++) {
        A.addEntry(i, i, diag_entry);
    }

    for (unsigned int i = 0; i < N - 1; i++) {
        A.addEntry(i, i + 1, upper_entry);
    }

    for (unsigned int i = 1; i < N; i++) {
        A.addEntry(i, i - 1, lower_entry);
    }


    //set up vector f
    Vector f(N, 0);
    f[N - 1] = -1 * upper_entry;
    f_ = f;

    A_ = A;


}


SparseMatrix FiniteDifferenceSolver::getMatrix() {
    return (A_);
}

Vector FiniteDifferenceSolver::getVector() {
    return (f_);
}

Vector FiniteDifferenceSolver::solver() {

    //set up initial guess for Gauss Seidel algorithm
    Vector x_0(gridPoints_N_, 0);

    //call Gauss seidel to invert A_ against f_

    return (A_.GaussSeidel(f_, x_0, 1e-6, 100000));


}

double FiniteDifferenceSolver::analyticSoln(double x, double alpha, double beta) {

    double c, ans;
    c = 1 / (-1 + exp(beta / alpha));

    ans = -1 * c + c * exp((beta / alpha) * x);

    return ans;


}

Vector FiniteDifferenceSolver::exactSolution() {

    //vector of the exact solution at the gridpoints

    double h = double(1) / double(gridPoints_N_ + 1);
    double x;

    Vector u(gridPoints_N_);

    for (unsigned int i = 0; i < u.size(); i++) {
        //u[0]= u_{1} = u(x_{1}) = u(1*h)
        x = (i + 1) * h;
        u[i] = analyticSoln(x, alpha_, beta_);
        cout << u[i] << endl;
    }

    return u;

}

double FiniteDifferenceSolver::vectorInfinityNorm(Vector vec) {

    vector<double>::iterator vecIterator;
    vecIterator = vec.begin();

    double temp_max = fabs((*vecIterator));

    for (vecIterator; vecIterator != vec.end(); vecIterator++) {

        if (fabs(*vecIterator) >= fabs(temp_max)) {
            temp_max = fabs(*vecIterator);
        }

    }

    return temp_max;

}

double FiniteDifferenceSolver::error(Vector x, Vector y) {

    Vector z(x.size());

    for (unsigned int i = 0; i < x.size(); i++) {
        z[i] = x[i] - y[i];
    }

    return vectorInfinityNorm(z);

}


Vector FiniteDifferenceSolver::getGridPoints(){

    Vector gridPoints(gridPoints_N_);

    for(unsigned int i=0; i< gridPoints.size(); i++)
    {

    }

}
