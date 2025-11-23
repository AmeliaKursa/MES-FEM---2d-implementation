import MESInputAndCalculations
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

        hbc = MESInputAndCalculations.HbcGenerator(
            element=element,
            grid=siatka,
            data=data1,
            el_univ=el_univ
        )

        hbc.compute()

        print("\nMacierz Hbc lokalna:")
        for r in range(4):
            print("  ".join(f"{element.Hbc[r][c]:10.6f}" for c in range(4)))

    print("\n==============================")
    print("AGREGACJA GLOBALNA")
    print("==============================")

    ukl = MESInputAndCalculations.UkladRownan(siatka)
    ukl.agregacja2d()

    print("\nMacierz H globalna:")
    for r in range(siatka.numNode):
        print("  ".join(f"{ukl.hGlobalna[r][c]:10.6f}" for c in range(siatka.numNode)))


if __name__ == "__main__":
    main()
