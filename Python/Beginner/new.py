from sympy import symbols

x, y = symbols('x,y')
expr = x + 2*y
print expr + 1
print expr - x

from sympy import expand, factor
expanded_expr = expand(x*expr)
print expanded_expr
print factor(expanded_expr)
