from sympy import init_printing, Rational

from classes.Bezier import Bezier

init_printing()

control_points = [[16, 0], [32, 32], [-32, 32], [0, 16], [-16, 0]]

bezier = Bezier(control_points)

bezier.raise_power()
