from sympy import *

from functions import (
    getLineAt,
    getPlane,
    getVectorDerivative,
    makeEquationChain,
    prettyPrintEquation,
    prettyPrintLine,
    prettyPrintPlane,
    prettyPrintVector,
    scalarProduct,
    vectorLengthMatrix,
    vectorProductMatrix,
)
from symbols import DELTA

init_printing()

# INPUT DATA HERE
u = Symbol("u")
v = Symbol("v")
a = Symbol("a", constant=True)
b = Symbol("b", constant=True)

rX = u * cos(v)
rY = u * sin(v)
rZ = b - a / u
c1U = v ** 2 + 3
c2U = u

# You may need to calculate these by hand if point is given as coords - e.g. (1,3,4)
atPointU = 1
atPointV = 2

# Calculate these by hand!!
cU = 1
cV = 2

rVector = [rX, rY, rZ]

r1d = rX.diff(u)

ruVectorD = getVectorDerivative(rVector, u)
rvVectorD = getVectorDerivative(rVector, v)

prettyPrintVector("-> r_u", ruVectorD)
prettyPrintVector("-> r_v", rvVectorD)

prod = vectorProductMatrix(ruVectorD, rvVectorD)
len = vectorLengthMatrix(prod)

prettyPrintVector("-> r_u x -> r_v", prod)
prettyPrintEquation("|-> r_u x -> r_v|", len)

nVector = prod / len
nVectorWithoutDenominator = prod

prettyPrintVector("-> N", nVector)

tangentLine = getLineAt(rVector, nVectorWithoutDenominator)

prettyPrintLine("N", tangentLine)

tangentPlane = getPlane(rVector, nVectorWithoutDenominator)

prettyPrintPlane("TpS", tangentPlane)

# TODO: do all subs calls in one call as an array
prettyPrintVector(
    f"-> N (M(u={atPointU}, v={atPointV}))", nVector.subs(u, atPointU).subs(v, atPointV)
)

prettyPrintLine(
    f"N(M(u={atPointU}, v={atPointV}))",
    [x.subs(u, atPointU).subs(v, atPointV) for x in tangentLine],
)

prettyPrintPlane(
    f"TpS(M(u={atPointU}, v={atPointV}))",
    tangentPlane.subs(u, atPointU).subs(v, atPointV),
)

g11 = Symbol("g_11")
g22 = Symbol("g_22")
g12 = Symbol("g_12")
g = Symbol("g")
du = Symbol("d_u")
dv = Symbol("d_v")


firstFormFormula = g11 * (du ** 2) + 2 * g12 * du * dv + g22 * (dv ** 2)

g11Value = scalarProduct(ruVectorD, ruVectorD)
g22Value = scalarProduct(rvVectorD, rvVectorD)
g12Value = scalarProduct(ruVectorD, rvVectorD)

gValue = g11Value * g22Value - g12Value ** 2

prettyPrintEquation(str(g11), g11Value)
prettyPrintEquation(str(g22), g22Value)
prettyPrintEquation(str(g12), g12Value)
prettyPrintEquation(str(g), gValue)

prettyPrintEquation(
    "I(du,dv)",
    makeEquationChain(
        [
            firstFormFormula,
            firstFormFormula.subs(g11, g11Value)
            .subs(g22, g22Value)
            .subs(g12, g12Value),
        ]
    ),
)

h11 = Symbol("h_11")
h22 = Symbol("h_22")
h12 = Symbol("h_12")
h = Symbol("h")


ruuVectorD = getVectorDerivative(rVector, u, 2)
rvvVectorD = getVectorDerivative(rVector, v, 2)
ruvVectorD = getVectorDerivative(getVectorDerivative(rVector, u), v)

prettyPrintVector("-> r_uu", ruuVectorD)
prettyPrintVector("-> r_vv", rvvVectorD)
prettyPrintVector("-> r_uv", ruvVectorD)

h11Value = scalarProduct(nVector, ruuVectorD)
h22Value = scalarProduct(nVector, rvvVectorD)
h12Value = scalarProduct(nVector, ruvVectorD)

hValue = h11Value * h22Value - h12Value ** 2

prettyPrintEquation(str(h11), h11Value)
prettyPrintEquation(str(h22), h22Value)
prettyPrintEquation(str(h12), h12Value)
prettyPrintEquation(str(h), hValue)

secondFormFormula = h11 * (du ** 2) + 2 * h12 * du * dv + h22 * dv ** 2

prettyPrintEquation(
    "II(du,dv)",
    makeEquationChain(
        [
            secondFormFormula,
            secondFormFormula.subs(h11, h11Value)
            .subs(h22, h22Value)
            .subs(h12, h12Value),
        ]
    ),
)

prettyPrintEquation(
    "II(du,dv) / I(du,dv)",
    (secondFormFormula / firstFormFormula).subs(
        [
            (du, 0.5),
            (dv, 1),
            (g11, g11Value),
            (g22, g22Value),
            (g12, g12Value),
            (h11, h11Value),
            (h22, h22Value),
            (h12, h12Value),
        ]
    ),
)

deltaU = Symbol(f"{DELTA}u")
deltaV = Symbol(f"{DELTA}v")

print("3.")

prettyPrintEquation("g_11^p", g11Value.subs(u, cU).subs(v, cV))
prettyPrintEquation("g_22^p", g22Value.subs(u, cU).subs(v, cV))
prettyPrintEquation("g_12^p", g12Value.subs(u, cU).subs(v, cV))

print("4.")

duValue = c1U.diff(v)
dvValue = u.diff(u)

prettyPrintEquation(
    "C_1",
    makeEquationChain([Symbol(f"u={c1U}"), Symbol(f"du=d({c1U})"), duValue, dvValue]),
)

print(f"C_1: (du, dv) :: ({duValue}, {dvValue})")
print()
print("5.")

deltaUValue = c2U.diff(v)
deltaVValue = u.diff(u)
prettyPrintEquation(
    "C_2",
    makeEquationChain(
        [
            Symbol(f"u={c2U}"),
            Symbol(f"{DELTA}u={DELTA}({c2U})"),
            deltaUValue,
            deltaVValue,
        ]
    ),
)
print(f"C_2: ({DELTA}u, {DELTA}v) :: ({deltaUValue}, {deltaVValue})")
print()

print("6.")
cosFormula = (
    g11 * du * deltaU + g12 * (du * deltaV + dv * deltaU) + g22 * dv * deltaV
) / (
    (sqrt(g11 * (du ** 2) + 2 * g12 * du * dv + g22 * (dv ** 2)))
    * (sqrt(g11 * (deltaU ** 2) + 2 * g12 * deltaU * deltaV + g22 * (deltaV ** 2)))
)

prettyPrintEquation("cos(C_1,C_2)", cosFormula)
prettyPrintEquation(
    "cos(C_1,C_2)",
    cosFormula.subs(g11, g11Value)
    .subs(g22, g22Value)
    .subs(g12, g12Value)
    .subs(du, duValue)
    .subs(dv, dvValue)
    .subs(deltaU, deltaUValue)
    .subs(deltaV, deltaVValue)
    .subs(u, cU)
    .subs(v, cV),
)

K = h / g

prettyPrintEquation("K", makeEquationChain([K, K.subs(g, gValue).subs(h, hValue)]))

H = (g11 * h22 - 2 * g12 * h12 + g22 * h11) / (2 * g)
prettyPrintEquation(
    "H",
    makeEquationChain(
        [
            H,
            simplify(
                H.subs(g, gValue)
                .subs(g11, g11Value)
                .subs(g12, g12Value)
                .subs(g22, g22Value)
                .subs(h11, h11Value)
                .subs(h12, h12Value)
                .subs(h22, h22Value)
            ),
        ]
    ),
)

# if hValue < 0:
#     print("2 asymptotic lines")
# elif hValue == 0:
#     print("1 asymptotic line")
# else:
#     print("No asymptotic lines exist")
