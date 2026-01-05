import MESInputAndCalculations
import numpy as np
import math


def main():
    siatka = MESInputAndCalculations.Grid()
    siatka.readGridData("Test1_4_4.txt")
    siatka.printGrid()

    data1 = MESInputAndCalculations.GlobalData()
    data1.readGlobalVariables("Test1_4_4.txt")
    data1.printGlobalVariables()

    el_univ = MESInputAndCalculations.ElementUniwersalny(npc=2)
    el_univ.pochodneKsi()
    el_univ.pochodneEta()

    for e_idx, element in enumerate(siatka.elementsArray):
        print(f"\n===== ELEMENT {e_idx + 1} =====")

        jakobian = MESInputAndCalculations.Jakobian(
            npc=2,
            element=element,
            grid=siatka,
            el_univ=el_univ,
            data=data1
        )
        jakobian.obliczanieWartosci()

        print("\nMacierz H lokalna:")
        for r in range(4):
            print("  ".join(f"{jakobian.H[r][c]:10.6f}" for c in range(4)))

        wbrzegowe = MESInputAndCalculations.HbcAndPGenerator(
            element=element,
            grid=siatka,
            data=data1,
            el_univ=el_univ
        )

        wbrzegowe.computeHbc()
        wbrzegowe.computeP()

        print("\nMacierz Hbc lokalna:")
        for r in range(4):
            print("  ".join(f"{element.Hbc[r][c]:10.6f}" for c in range(4)))

        print("\nWektor P lokalny:")
        for r in range(4):
            print(f"{element.P[r]:10.3f}")

        print("\nMacierz C lokalna:")
        for r in range(4):
            print("  ".join(f"{element.C[r][c]:10.6f}" for c in range(4)))

    print("\n==============================")
    print("AGREGACJA GLOBALNA")
    print("==============================")

    ukl = MESInputAndCalculations.UkladRownan(siatka)
    ukl.agregacja2dH()
    ukl.agregacja2dC()
    ukl.agregacja2dP()

    print("\nMacierz H globalna:")
    for r in range(siatka.numNode):
        print("  ".join(f"{ukl.hGlobalna[r][c]:10.6f}" for c in range(siatka.numNode)))

    print("\nMacierz C globalna:")
    for r in range(siatka.numNode):
        print("  ".join(f"{ukl.cGlobalna[r][c]:10.6f}" for c in range(siatka.numNode)))

    print("\nWektor P globalny:")
    for r in range(siatka.numNode):
        print(f"{ukl.pGlobalne[r]:10.3f}")

    print("\n=========ROZWIAZYWANIE UKLADU ROWNAN (NIEUSTALONE)=========")

    H = np.array(ukl.hGlobalna, dtype=float)
    C = np.array(ukl.cGlobalna, dtype=float)
    P = np.array(ukl.pGlobalne, dtype=float)

    T = np.array(
        [data1.InitialTemp for _ in range(siatka.numNode)],
        dtype=float
    )

    dt = data1.SimulationStepTime
    n_steps = int(data1.SimulationTime / dt)

    for step in range(1, n_steps + 1):
        time = step * dt

        A = H + C / dt
        b = P + (C / dt) @ T

        T = np.linalg.solve(A, b)

        print(
            f"t = {time:6.1f}  "
            f"Tmin = {T.min():10.5f}  "
            f"Tmax = {T.max():10.5f}"
        )


if __name__ == "__main__":
    main()
