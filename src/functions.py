from sympy import *


def normalizeInfinity(expr):
    divisionByZero = 1 / Symbol("O")
    # zoo is complex infinity, which is the result of 1 / 0.
    # We replace it with 1 / O (this is the 15th letter of the alphabet, not
    # zero) so that it looks better visually.
    return expr.subs(zoo, divisionByZero)


def normalizeArray(array):
    return [normalizeInfinity(x) for x in array]


def vectorLength(a):
    len = sqrt(a[0] ** 2 + a[1] ** 2 + a[2] ** 2)
    return simplify(len)


def vectorProduct(a, b):
    matrix1 = Matrix([[a[1], a[2]], [b[1], b[2]]])
    matrix2 = Matrix([[a[0], a[2]], [b[0], b[2]]])
    matrix3 = Matrix([[a[0], a[1]], [b[0], b[1]]])

    det1 = matrix1.det()
    det2 = matrix2.det()
    det3 = matrix3.det()

    simplified = [simplify(x) for x in [det1, -det2, det3]]
    return simplified


def scalarProduct(vectorA, vectorB):
    return simplify(
        vectorA[0] * vectorB[0] + vectorA[1] * vectorB[1] + vectorA[2] * vectorB[2]
    )


def prettyPrintEquation(label, equation):
    normalized = normalizeInfinity(equation)
    pprint(Eq(Symbol(label), normalized, evaluate=False))
    print()


def prettyPrintVector(label, vector):
    normalized = normalizeArray(vector)
    vector = Matrix(normalized).transpose()
    pprint(Eq(Symbol(label), vector, evaluate=False))
    print()


def verifyVector(vector, label):
    # After normalization maybe we can remove this?
    for i in range(len(vector)):
        if vector[i] == 0:
            print(
                f"WARNING: during calculation of {label}, vector[{i}] is 0 so there may be errors"
            )


def getLineAt(point, vector, label=None):
    x = Symbol("x")
    y = Symbol("y")
    z = Symbol("z")

    verifyVector(vector, label)

    equation1 = (x - point[0]) / vector[0]
    equation2 = (y - point[1]) / vector[1]
    equation3 = (z - point[2]) / vector[2]

    lineVector = [simplify(x) for x in [equation1, equation2, equation3]]

    return lineVector


def prettyPrintLine(label, line):
    LAMBDA = "\u03BB"
    normalized = normalizeArray(line)
    pprint(label)
    pprint(
        Eq(
            Eq(
                Eq(normalized[0], normalized[1], evaluate=False),
                normalized[2],
                evaluate=False,
            ),
            Symbol(LAMBDA),
            evaluate=False,
        )
    )
    print()


def prettyPrintPlane(label, plane):
    pprint(label)
    pprint(Eq(plane, 0))
    print()


def getIntersectionLine(label, line, replacement, point):
    replacedLine = []
    for i in range(len(line)):
        replacedLine.append(line[i].subs(replacement, point))

    prettyPrintLine(label, replacedLine)


def getPlane(lineVector, point):
    x = Symbol("x")
    y = Symbol("y")
    z = Symbol("z")

    plane = (
        point[0] * (x - lineVector[0])
        + point[1] * (y - lineVector[1])
        + point[2] * (z - lineVector[2])
    )

    return plane


def getIntersectionPlane(label, plane, replacement, point):
    replacedPlane = plane.subs(replacement, point)
    prettyPrintPlane(label, replacedPlane)


def getVectorDerivative(vector, variable, derivative=1):
    return [x.diff(variable, derivative) for x in vector]


def getKappa(vector, variable):
    vectorD = getVectorDerivative(vector, variable)
    vectorDD = getVectorDerivative(vector, variable, 2)

    kappa = vectorLength(vectorProduct(vectorD, vectorDD)) / (
        vectorLength(vectorD) ** 3
    )

    return kappa


def getTau(vector, variable):
    vectorD = getVectorDerivative(vector, variable)
    vectorDD = getVectorDerivative(vector, variable, 2)
    vectorDDD = getVectorDerivative(vector, variable, 3)

    tau = scalarProduct(vectorProduct(vectorD, vectorDD), vectorDDD) / (
        vectorLength(vectorProduct(vectorD, vectorDD)) ** 2
    )
    return tau
