from sympy import *

from functions import prettyPrintEquation

init_printing()


controlPoints = [[3, 0], [0, 3], [-3, 0], [0, -3]]

uValue = Rational("3/10")

curvePower = len(controlPoints) - 1

inverseU = 1 - uValue

A = Symbol("A")
B = Symbol("B")

C = inverseU * A + uValue * B

points = [Matrix(cp).transpose() for cp in controlPoints]
for i in range(len(points)):
    prettyPrintEquation(f"P{i}", points[i])

for i in range(1, curvePower + 1):
    newPoints = []
    for j in range(curvePower - i + 1):
        newPoint = C.subs(A, points[j]).subs(B, points[j + 1])
        newPoints.append(newPoint)
        prettyPrintEquation(f"P{i}{j}", newPoint)

    points = newPoints
    pass
