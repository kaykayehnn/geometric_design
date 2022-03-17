from sympy import Symbol, nan, Integer
from functions import (
    convertControlPointsToVectors,
    findPositionInBSpline,
    incrementLetter,
    incrementLetterPrimeness,
    makeEquationChain,
    prettifySymbol,
    prettyPrintEquation,
    printControlPointsTable,
    printKnots,
    trimTrailingComma,
    verifyBSpline,
)
from symbols import IDENTITY, EPSILON


class BSpline:
    def __init__(self, knots):
        self.knots = knots

    def calculate(self, u_value, target_power):
        curvePower = len(self.knots) - 1

        u = Symbol("u")
        i = Symbol("i")

        for i in range(curvePower):
            ui = prettifySymbol(f"u_{i}")
            ui1 = prettifySymbol(f"u_{i+1}")
            n = prettifySymbol(f"N_{i}0")
            closingBracket = "]" if self.knots[i + 1] == 1 else ")"
            prefix = ""
            if self.knots[i] < u_value and self.knots[i + 1] > u_value:
                prefix = f"{u_value} {EPSILON} "
            print(
                f"{prefix}[ {ui} ,  {ui1} ) {IDENTITY} [ {self.knots[i]} ; {self.knots[i+1]} {closingBracket} -> {n}"
            )

        print()

        # Calculate first row
        nValues = []
        for i in range(curvePower):
            value = Integer(
                1 if self.knots[i] < u_value and self.knots[i + 1] > u_value else 0
            )
            nValues.append(value)
            prettyPrintEquation(f"{prettifySymbol(f'N_{i}0')} ({u_value})", value)

        for p in range(1, target_power + 1):
            newNValues = []
            for i in range(curvePower - p):
                nip = prettifySymbol(f"N{i}{p}") + "(u)"
                ui = Symbol(f"u{i}")
                uip = Symbol(f"u{i+p}")
                uip1 = Symbol(f"u{i+p+1}")
                ui1 = Symbol(f"u{i+1}")
                nipM1 = Symbol(f"N{i},{p-1}")
                ni1pM1 = Symbol(f"N{i+1},{p-1}")
                formula = (u - ui) / (uip - ui) * nipM1 + (uip1 - u) / (
                    uip1 - ui1
                ) * ni1pM1

                symbolicFormula = formula.subs(u, u_value)

                tempFormula = (
                    formula.subs(ui, self.knots[i])
                    .subs(uip, self.knots[i + p])
                    .subs(uip1, self.knots[i + p + 1])
                    .subs(ui1, self.knots[i + 1])
                    .subs(u, u_value)
                )

                firstStageFormula = tempFormula.subs(nipM1, Symbol(f"N{i}{p-1}")).subs(
                    ni1pM1, Symbol(f"N{i+1}{p-1}")
                )
                result = tempFormula.subs(nipM1, nValues[i]).subs(
                    ni1pM1, nValues[i + 1]
                )

                output = makeEquationChain([symbolicFormula, firstStageFormula, result])

                prettyPrintEquation(nip, output)
                newNValues.append(result if result != nan else 0)

            nValues = newNValues

    def insert(self, control_points, points_to_insert):
        cp = Symbol("i")
        t = Symbol("t")

        START_KNOT_LETTER = "u"
        START_CP_LETTER = "P"

        current_knots = self.knots
        controlPointVectors = convertControlPointsToVectors(control_points)
        curvePower = len(current_knots) - len(control_points) - 1
        knotLetter = START_KNOT_LETTER
        # controlPointLetter stores the P, P', P''... letters
        controlPointLetter = START_CP_LETTER
        # tempControlPointLetter stores the P, Q, R... letters
        tempControlPointLetter = START_CP_LETTER

        for i in range(len(points_to_insert)):
            tValue = points_to_insert[i]
            aiSuffix = "" if i == 0 else str(i + 1)
            ai = Symbol("a_i" + aiSuffix)

            verifyBSpline(current_knots)

            printKnots(current_knots, knotLetter)

            index = findPositionInBSpline(current_knots, tValue)
            affectedPointsOutput = f"t {EPSILON} [ {prettifySymbol(knotLetter + str(index))}, {prettifySymbol(knotLetter + str(index+1))} ) -> "

            # It is important to note that these values correspond to the P indexes.
            affectedPointsStart = index - curvePower
            affectedPointsEnd = index

            separator = ", "
            for cp in range(affectedPointsEnd, affectedPointsStart - 1, -1):
                affectedPointsOutput += (
                    prettifySymbol(f"{controlPointLetter}_{cp}") + separator
                )

            affectedPointsOutput = trimTrailingComma(
                affectedPointsOutput, len(separator)
            )
            print(affectedPointsOutput)
            print()

            affectedPoints = []
            newTempCPLetter = incrementLetter(tempControlPointLetter)
            qPointIndexStart = affectedPointsStart + 1
            qPointIndexEnd = affectedPointsEnd + 1

            newControlPoints = controlPointVectors[:qPointIndexStart]

            aiValues = []
            for cp in range(qPointIndexStart, qPointIndexEnd):
                ui = Symbol(f"u_{cp}")
                uip = Symbol(f"u_{cp+curvePower}")
                aiFormula = (t - ui) / (uip - ui)
                aiValue = (
                    aiFormula.subs(t, tValue)
                    .subs(ui, current_knots[cp])
                    .subs(uip, current_knots[cp + curvePower])
                )
                aiOutput = makeEquationChain(
                    [
                        aiFormula,
                        aiValue,
                    ]
                )

                prettyPrintEquation(f"a_{cp}{aiSuffix}", aiOutput)
                aiValues.append(aiValue)

            for cp in range(qPointIndexStart, qPointIndexEnd):
                piM1 = Symbol(f"{controlPointLetter}_{cp-1}")
                pi = Symbol(f"{controlPointLetter}_{cp}")
                qiFormula = (1 - ai) * piM1 + ai * pi

                aiValue = aiValues[cp - affectedPointsStart - 1]

                qiValue = (
                    qiFormula.subs(ai, aiValue)
                    .subs(piM1, controlPointVectors[cp - 1])
                    .subs(pi, controlPointVectors[cp])
                )
                qiOutput = makeEquationChain([qiFormula, qiValue.transpose()])

                prettyPrintEquation(f"{newTempCPLetter}_{cp}", qiOutput)
                newControlPoints.append(qiValue)

            newControlPoints.extend(controlPointVectors[affectedPointsEnd:])

            newControlPointVectors = convertControlPointsToVectors(newControlPoints)
            newKnotLetter = incrementLetter(knotLetter)
            newControlPointLetter = incrementLetterPrimeness(controlPointLetter)
            newKnots = current_knots.copy()
            newKnots.insert(index + 1, tValue)

            printKnots(newKnots, newKnotLetter)

            printControlPointsTable(
                controlPointLetter,
                len(newControlPoints),
                newTempCPLetter,
                qPointIndexStart,
                qPointIndexEnd,
            )

            knotLetter = newKnotLetter
            controlPointLetter = newControlPointLetter
            tempControlPointLetter = newTempCPLetter
            current_knots = newKnots
            controlPointVectors = newControlPointVectors

            if i != len(points_to_insert) - 1:
                # Print separator before starting next insert
                print("------------------------------------------")
                print()
