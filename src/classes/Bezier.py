from sympy import (
    Symbol,
    factorial,
    pretty,
    Eq,
    simplify,
    Matrix,
    UnevaluatedExpr,
    Rational,
)
from functions import (
    convertControlPointsToVectors,
    displayCurve,
    displayCurveRaisePowerFormula,
    makeEquationChain,
    prettifySymbol,
    prettyPrintEquation,
)


class Bezier:
    def __init__(self, control_points):
        self.control_points = control_points

    def calculate(self, u_value):
        curvePower = len(self.control_points) - 1

        u = Symbol("u")
        n = Symbol("n")
        i = Symbol("i")

        # B(n,i)
        B = (
            factorial(n)
            / (factorial(i) * factorial(n - i))
            * (u ** i)
            * (1 - u) ** (n - i)
        )

        equations = []

        outputs = []
        # Only used for error reporting
        sum = 0

        for j in range(curvePower + 1):
            equation = B.subs(n, curvePower).subs(i, j)
            value = equation.subs(u, u_value)

            if value < 0 or value > 1:
                raise Exception(
                    f"SOMETHING WENT WRONG: value must be [0, 1] but is {value}"
                )

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
            xEquation += equation * self.control_points[j][0]
            yEquation += equation * self.control_points[j][1]

        xEquation = simplify(xEquation)
        yEquation = simplify(yEquation)

        prettyPrintEquation("C(u)", C)
        prettyPrintEquation("x(u)", xEquation)
        prettyPrintEquation("y(u)", yEquation)
        prettyPrintEquation("C(u)", Matrix([xEquation, yEquation]).transpose())

        result = [x.subs(u, u_value) for x in [xEquation, yEquation]]

        prettyPrintEquation(f"C(u0={u_value})", Matrix(result).transpose())

    def casteljau(self, u_value):
        curvePower = len(self.control_points) - 1

        inverseU = 1 - u_value

        u = Symbol("u")

        points = convertControlPointsToVectors(self.control_points)
        for i in range(len(points)):
            prettyPrintEquation(f"P{i}", points[i])

        allPoints = [points]

        for i in range(1, curvePower + 1):
            newPoints = []
            for j in range(curvePower - i + 1):
                A = Symbol(f"P{i-1}{j}")
                B = Symbol(f"P{i-1}{j+1}")
                equation = (1 - u) * A + u * B
                replacedEquation = (
                    equation.subs(u, UnevaluatedExpr(u_value))
                    .subs(A, UnevaluatedExpr(points[j]))
                    .subs(B, UnevaluatedExpr(points[j + 1]))
                )
                result = replacedEquation.doit()

                finalFormula = Eq(
                    equation,
                    Eq(replacedEquation, result, evaluate=False),
                    evaluate=False,
                )
                newPoints.append(result)
                prettyPrintEquation(f"P{i}{j}", finalFormula)

            points = newPoints
            allPoints.append(newPoints)

        # Calculate first and second derivative
        secondToLastLine = allPoints[curvePower - 1]
        C1 = curvePower * (secondToLastLine[1] - secondToLastLine[0])

        prettyPrintEquation(f"C.({u_value})", C1)

        thirdToLastLine = allPoints[curvePower - 2]
        C2 = (
            curvePower
            * (curvePower - 1)
            * (thirdToLastLine[0] - 2 * thirdToLastLine[1] + thirdToLastLine[2])
        )

        prettyPrintEquation(f"C..({u_value})", C2)

    def raise_power(self):
        u = Symbol("u")
        i = Symbol("i")

        curvePower = len(self.control_points) - 1

        oldCurve = "C(u)"
        newCurve = "D(u)"
        oldPoint = "P"
        newPoint = "Q"
        newCurvePower = curvePower + 1
        newTotalPoints = newCurvePower + 1

        # Display the current curve
        displayCurve(oldCurve, self.control_points)
        # Display the new curve formula
        displayCurveRaisePowerFormula(newCurve, newPoint, curvePower)
        # Display the start and end points, which are the same
        print(
            f'{prettifySymbol(newPoint + "0")}={prettifySymbol(oldPoint + "0")} => {prettifySymbol(newPoint + "0")}{pretty(self.control_points[0])}'
        )
        print(
            f"{prettifySymbol(newPoint + str(newCurvePower))}={prettifySymbol(oldPoint + str(curvePower))} => {prettifySymbol(newPoint + str(newCurvePower))}{pretty(self.control_points[curvePower])}"
        )

        controlPointVectors = convertControlPointsToVectors(self.control_points)
        for i in range(1, newTotalPoints - 1):
            A = Symbol(oldPoint + str(i - 1))
            B = Symbol(oldPoint + str(i))
            Q = (
                UnevaluatedExpr(Rational(i, newCurvePower)) * A
                + (1 - UnevaluatedExpr(Rational(i, newCurvePower))) * B
            )
            Qreplaced = Q.subs(A, UnevaluatedExpr(controlPointVectors[i - 1])).subs(
                B, UnevaluatedExpr(controlPointVectors[i])
            )
            result = Qreplaced.doit()

            equation = makeEquationChain([Q, Qreplaced, result])

            prettyPrintEquation(newPoint + str(i), equation)
