#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>

#ifndef C1_ASSIGNMENT3_NEWTONRAPH_H
#define C1_ASSIGNMENT3_NEWTONRAPH_H

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>
#include "models.hh"

template<class Model>
double
newton(double time, const double &y, double h, const Model &model, double a, double sum, double c_0) {
    double const tol = 1e-6;
    int const maxIter = 100000;
    for (int i = 0; i < maxIter; i++) {
        double x = c_0;
//        double fx = model.f(t, x) - ((x - y - h * sum) / h * a);
//        double fx1 = model.df(t, x) - (1 / h * a);
        c_0 = x - (double) (model.f(time, y + h * sum + h * a * c_0) - c_0) /
                  (h * a * model.df(time, y + h * sum + h * a * c_0) - 1);


        if (fabs(c_0 - x) < tol) {
            return c_0;
        }
    }

    std::cout << "max iterations hit" << std::endl;
    return c_0;
}


#endif //C1_ASSIGNMENT3_NEWTONRAPH_H
