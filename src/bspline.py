from sympy import *

from classes.BSpline import BSpline

init_printing()

knots = [
    0,
    0,
    0,
    Rational(2, 5),
    Rational(1, 2),
    Rational(3, 5),
    1,
    1,
    1,
]
u_value = Rational(7, 10)
target_power = 2

bspline = BSpline(knots)

bspline.calculate(u_value, target_power)
