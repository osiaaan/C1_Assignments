from sympy import *
from sympy import init_printing
init_printing() # doctest: +SKIP
x, t, z, nu = symbols('x t z nu')

init_printing(use_unicode=True)

init_printing()
print diff(sin(pi*x/(0.25 + x*y)),x)
print integrate(exp(x)*sin(x) + exp(x)*cos(x), x)
print integrate(sqrt(1/x),x)
print integrate(sin(x**2), (x,-oo,oo))
print simplify(sin(x)**2+cos(x)**2)
print factor(x**2 + 2*x + 1)
