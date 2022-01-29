from sympy import *
from functions import (
    getIntersectionLine,
    getIntersectionPlane,
    getKappa,
    getLineAt,
    getPlane,
    getTau,
    getVectorDerivative,
    prettyPrintLine,
    prettyPrintPlane,
    prettyPrintVector,
    vectorLength,
    vectorProduct,
)

init_printing()

x = Symbol("x")
u = Symbol("u")

a = Symbol("a", real=True, positive=True)
b = Symbol("b", real=True, positive=True)
c = Symbol("c")

rX = a * cos(u)
rY = a * sin(u)
rZ = b * u

rVector = [rX, rY, rZ]

ruVector = getVectorDerivative(rVector, u)
[rXu, rYu, rZu] = ruVector

prettyPrintVector("r. =", ruVector)

r1Len = vectorLength(ruVector)

pprint("|r.|=")
pprint(r1Len)

ruuVector = getVectorDerivative(rVector, u, 2)
[rXuu, rYuu, rZuu] = ruVector

prettyPrintVector("r.. =", ruuVector)

# Find vector product r. x r1..

productVector = vectorProduct(ruVector, ruuVector)

prettyPrintVector("r. x r.. =", productVector)

productVectorLen = vectorLength(productVector)

pprint("|r. x r..| =")
pprint(productVectorLen)

tVector = [x / r1Len for x in ruVector]

prettyPrintVector("t-> =", tVector)

bVectorX = simplify(productVector[0] / productVectorLen)
bVectorY = simplify(productVector[1] / productVectorLen)
bVectorZ = simplify(productVector[2] / productVectorLen)

bVector = [bVectorX, bVectorY, bVectorZ]

prettyPrintVector("b-> =", [bVectorX, bVectorY, bVectorZ])

mVector = vectorProduct(bVector, tVector)

prettyPrintVector("m -> =", mVector)

# Here we use only the numerator part of the vectors for simplification (all 3
# fractions share the same denominator so, we can remove it)
tLine = getLineAt(rVector, ruVector, "t")
bLine = getLineAt(rVector, productVector, "b")
mLine = getLineAt(rVector, mVector, "m")

prettyPrintLine("t : ", tLine)
prettyPrintLine("b : ", bLine)
prettyPrintLine("m : ", mLine)

MPoint = 0
if MPoint is not None:
    getIntersectionLine("t (M(u={MPoint}) :", tLine, u, MPoint)
    getIntersectionLine("b (M(u={MPoint}) :", bLine, u, MPoint)
    getIntersectionLine("m (M(u={MPoint}) :", mLine, u, MPoint)


NU = "\u03BD"
MU = "\u03BC"
PI = "\u03C0"

nuPlane = getPlane(rVector, ruVector)
muPlane = getPlane(rVector, productVector)
piPlane = getPlane(rVector, mVector)

prettyPrintPlane(f"{NU} : ", nuPlane)
prettyPrintPlane(f"{MU} : ", muPlane)
prettyPrintPlane(f"{PI} : ", piPlane)

MPointPlane = 0
if MPointPlane is not None:
    getIntersectionPlane(f"{NU} (u={MPointPlane}) : ", nuPlane, u, MPointPlane)
    getIntersectionPlane(f"{MU} (u={MPointPlane}) : ", muPlane, u, MPointPlane)
    getIntersectionPlane(f"{PI} (u={MPointPlane}) : ", piPlane, u, MPointPlane)


KAPPA = "\u00E6"
TAU = "\u03C4"

rXuuu = rX.diff(u, 3)
rYuuu = rY.diff(u, 3)
rZuuu = rZ.diff(u, 3)

ruuuVector = [rX, rY, rZ]

kappaEquation = getKappa(rVector, u)
tauEquation = getTau(rVector, u)

pprint(f"{KAPPA} = ")
pprint(kappaEquation)
print()

pprint(f"{TAU} = ")
pprint(tauEquation)
