from sympy import *
from functions import (
    getIntersectionLine,
    getIntersectionPlane,
    getKappa,
    getLineAt,
    getPlane,
    getTau,
    getVectorDerivative,
    prettyPrintEquation,
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

rX = u * sqrt(2)
rY = (u ** 2) / 2
rZ = ln(u)

rVector = [rX, rY, rZ]

ruVector = getVectorDerivative(rVector, u)
[rXu, rYu, rZu] = ruVector

prettyPrintVector("r.", ruVector)

r1Len = vectorLength(ruVector)

prettyPrintEquation("|r.|", r1Len)

ruuVector = getVectorDerivative(rVector, u, 2)
[rXuu, rYuu, rZuu] = ruVector

prettyPrintVector("r..", ruuVector)

# Find vector product r. x r1..

productVector = vectorProduct(ruVector, ruuVector)

prettyPrintVector("r. x r..", productVector)

productVectorLen = vectorLength(productVector)

prettyPrintEquation("|r. x r..|", productVectorLen)

tVector = [x / r1Len for x in ruVector]

prettyPrintVector("t->", tVector)

bVectorX = simplify(productVector[0] / productVectorLen)
bVectorY = simplify(productVector[1] / productVectorLen)
bVectorZ = simplify(productVector[2] / productVectorLen)

bVector = [bVectorX, bVectorY, bVectorZ]

# TODO: https://stackoverflow.com/a/56653481/6317168
prettyPrintVector("b->", [bVectorX, bVectorY, bVectorZ])

mVector = vectorProduct(bVector, tVector)

prettyPrintVector("m->", mVector)

# Sanity check that the of Frene's vectors are all 1
if (
    vectorLength(tVector) != 1
    or vectorLength(bVector) != 1
    or vectorLength(mVector) != 1
):
    raise Exception("INVARIANT: Frene's vector have to have a length of 1")

# Here we use only the numerator part of the vectors for simplification (all 3
# fractions share the same denominator so, we can remove it)
tLine = getLineAt(rVector, ruVector, "t")
bLine = getLineAt(rVector, productVector, "b")
mLine = getLineAt(rVector, mVector, "m")

prettyPrintLine("t : ", tLine)
prettyPrintLine("b : ", bLine)
prettyPrintLine("m : ", mLine)

MPoint = 1
if MPoint is not None:
    getIntersectionLine(f"t (M(u={MPoint}) :", tLine, u, MPoint)
    getIntersectionLine(f"b (M(u={MPoint}) :", bLine, u, MPoint)
    getIntersectionLine(f"m (M(u={MPoint}) :", mLine, u, MPoint)


NU = "\u03BD"
MU = "\u03BC"
PI = "\u03C0"

nuPlane = getPlane(rVector, ruVector)
muPlane = getPlane(rVector, productVector)
piPlane = getPlane(rVector, mVector)

prettyPrintPlane(f"{NU} : ", nuPlane)
prettyPrintPlane(f"{MU} : ", muPlane)
prettyPrintPlane(f"{PI} : ", piPlane)

MPointPlane = 1
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

prettyPrintEquation(KAPPA, kappaEquation)
prettyPrintEquation(TAU, tauEquation)
