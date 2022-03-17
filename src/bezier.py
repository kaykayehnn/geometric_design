from sympy import init_printing, Rational
from classes.Bezier import Bezier

init_printing()

control_points = [
    [-2, 0],
    [0, 1],
    [3, 0],
]

u_value = Rational(1, 4)

b = Bezier(control_points)

b.calculate(u_value)
