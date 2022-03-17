from sympy import init_printing, Rational
from classes.BSpline import BSpline

init_printing()

# INPUT DATA HERE
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
# fmt: off
control_points = [
  [2,2],
  [0,2],
  [0,0],
  [-2,0],
  [0,-2],
  [2,-2],
]
# fmt:on
points_to_insert = [
    Rational(7, 10),
    Rational(7, 10),
]

bspline = BSpline(knots)

bspline.insert(control_points, points_to_insert)
