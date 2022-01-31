from sympy import *

from functions import prettyPrintEquation

init_printing()

controlPoints = [[3, 0], [0, 3], [-3, 0], [0, -3]]

uValue = Rational("3/10")

curvePower = len(controlPoints) - 1

u = Symbol("u")
n = Symbol("n")
i = Symbol("i")

# B(n,i)
B = factorial(n) / (factorial(i) * factorial(n - i)) * (u ** i) * (1 - u) ** (n - i)

equations = []

outputs = []
# Only used for error reporting
sum = 0

for j in range(curvePower + 1):
    equation = B.subs(n, curvePower).subs(i, j)
    value = equation.subs(u, uValue)

    if value < 0 or value > 1:
        raise Exception(f"SOMETHING WENT WRONG: value must be [0, 1] but is {value}")

    index = f"B({curvePower},{j})"

    evaluated = Eq(
        Symbol(f"{index}"), Eq(equation, value, evaluate=False), evaluate=False
    )
    output = pretty(evaluated)

    sum += value
    equations.append(equation)
    outputs.append(output)

if sum != 1:
    raise Exception(f"SOMETHING WENT WRONG: sum must be 1 but is {sum}")


for output in outputs:
    print(output)

C = 0
xEquation = 0
yEquation = 0
for j in range(len(equations)):
    equation = equations[j]
    C += equation * Symbol(f"P{j}")
    xEquation += equation * controlPoints[j][0]
    yEquation += equation * controlPoints[j][1]

xEquation = simplify(xEquation)
yEquation = simplify(yEquation)

prettyPrintEquation("C(u)", C)
prettyPrintEquation("C(u)", Matrix([xEquation, yEquation]).transpose())
prettyPrintEquation("x(u)", xEquation)
prettyPrintEquation("y(u)", yEquation)

result = [x.subs(u, uValue) for x in [xEquation, yEquation]]

prettyPrintEquation(f"C(u0={uValue})", Matrix(result).transpose())
