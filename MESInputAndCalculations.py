import math

from MESgauss import Gauss


def dNx_po_dksi(eta, n):
    if n == 0:
        return -0.25 * (1 - eta)
    elif n == 1:
        return 0.25 * (1 - eta)
    elif n == 2:
        return 0.25 * (1 + eta)
    elif n == 3:
        return -0.25 * (1 + eta)
    else:
        raise ValueError("N jest w przedziale <0,3>")


def dNx_po_deta(ksi, n):
    if n == 0:
        return -0.25 * (1 - ksi)
    elif n == 1:
        return -0.25 * (1 + ksi)
    elif n == 2:
        return 0.25 * (1 + ksi)
    elif n == 3:
        return 0.25 * (1 - ksi)
    else:
        raise ValueError("N jest w przedziale <0,3>")


class Node:

    def __init__(self, id1, xo, yo):
        self.x = xo
        self.y = yo
        self.id=id1
        self.BC=False



class Element:
    def __init__(self, v1, v2, v3, v4):
        self.nodeArray = [v1, v2, v3, v4]
        self.jakobian = []
        self.H = []
        self.Hbc=[]
        self.surfaces=[Surface(), Surface(), Surface(), Surface()]

class Surface:
    def __init__(self):
        self.pc=[]
        self.N=[]


class Grid:
    def __init__(self):
        self.numNode = 0
        self.numElem = 0
        self.nodesArray = []
        self.elementsArray = []
        self.bcArray = []

    def readGridData(self, filename):
        with open(filename, 'r') as f:
            lines = f.readlines()
            startReadingNode = False
            startReadingElement = False
            startReadingBC = False
            for line in lines:
                line = line.strip()
                if not line:
                    continue

                if line.startswith('*Node'):
                    startReadingNode = True
                    startReadingElement = False
                    startReadingBC = False
                    continue

                if line.startswith('*Element'):
                    startReadingNode = False
                    startReadingElement = True
                    startReadingBC = False
                    continue

                if line.startswith('*BC'):
                    startReadingNode = False
                    startReadingElement = False
                    startReadingBC = True
                    continue

                if startReadingNode:
                    parts = line.split(',')
                    xy = [float(x) for x in parts[0:]]
                    self.nodesArray.append(Node(xy[0], xy[1], xy[2]))

                elif startReadingElement:
                    parts = line.split(',')
                    nodes = [int(x) for x in parts[1:]]
                    self.elementsArray.append(Element(nodes[0], nodes[1], nodes[2], nodes[3]))

                elif startReadingBC:
                    parts = line.split(',')
                    nodes = [int(x) for x in parts[0:]]
                    self.bcArray.append(nodes)
            for i in range(len(self.bcArray[0])):
                for j in range(len(self.nodesArray)):
                    if self.bcArray[0][i] == int(self.nodesArray[j].id):
                        self.nodesArray[j].BC = True

            self.numNode = len(self.nodesArray)
            self.numElem = len(self.elementsArray)

    def printGrid(self):
        print('Node Array:')
        for node in self.nodesArray:
            print(f'{int(node.id)}: {node.x}, {node.y} - czy bc: {node.BC}')
        print('Elements Array:')
        for i, element in enumerate(self.elementsArray):
            print(f'{i + 1}: {element.nodeArray}')
        print('BC Array:')
        for i, bc in enumerate(self.bcArray):
            print(f'{i + 1}: {bc}')



class GlobalData:
    def __init__(self):
        self.SimulationTime = 0
        self.SimulationStepTime = 0
        self.Conductivity = 0
        self.Alfa = 0
        self.Tot = 0
        self.InitialTemp = 0
        self.Density = 0
        self.SpecificHeat = 0
        self.nN = 0
        self.nE = 0

    def readGlobalVariables(self, fileName):
        with open(fileName, 'r+') as file:
            for line in file:
                words = line.strip().split()
                for i, word in enumerate(words):
                    if word == 'SimulationTime':
                        self.SimulationTime = float(words[i + 1])
                    if word == 'SimulationStepTime':
                        self.SimulationStepTime = float(words[i + 1])
                    if word == 'Conductivity':
                        self.Conductivity = float(words[i + 1])
                    if word == 'Alfa':
                        self.Alfa = float(words[i + 1])
                    if word == 'Tot':
                        self.Tot = float(words[i + 1])
                    if word == 'InitialTemp':
                        self.InitialTemp = float(words[i + 1])
                    if word == 'Density':
                        self.Density = float(words[i + 1])
                    if word == 'SpecificHeat':
                        self.SpecificHeat = float(words[i + 1])
                    if word == 'number' and words[i - 1] == 'Nodes':
                        self.nN = int(words[i + 1])
                    if word == 'number' and words[i - 1] == 'Elements':
                        self.nE = int(words[i + 1])

    def printGlobalVariables(self):
        print('SimulationTime:', self.SimulationTime)
        print('SimulationStepTime:', self.SimulationStepTime)
        print('Conductivity:', self.Conductivity)
        print('Alfa:', self.Alfa)
        print('Tot:', self.Tot)
        print('Initial Temp:', self.InitialTemp)
        print('Density:', self.Density)
        print('Specific Heat:', self.SpecificHeat)
        print('nN:', self.nN)
        print('nE:', self.nE)


class ElementUniwersalny:
    def __init__(self, npc):
        self.npc = npc
        self.npc2d = npc * npc
        self.dN_dKsi = [[0.0 for _ in range(self.npc2d)] for _ in range(4)]
        self.dN_dEta = [[0.0 for _ in range(self.npc2d)] for _ in range(4)]
        if npc == 1:
            self.pc = [0]
            self.w = [2]
        elif npc == 2:
            self.pc = [-1 / math.sqrt(3), 1 / math.sqrt(3)]
            self.w = [1, 1]
        elif npc == 3:
            self.pc = [-math.sqrt(3 / 5), 0, math.sqrt(3 / 5)]
            self.w = [5 / 9, 8 / 9, 5 / 9]
        elif npc == 4:
            a1 = math.sqrt(3 / 7 - 2 / 7 * math.sqrt(6 / 5))
            a2 = math.sqrt(3 / 7 + 2 / 7 * math.sqrt(6 / 5))
            w1 = (18 + math.sqrt(30)) / 36
            w2 = (18 - math.sqrt(30)) / 36
            self.pc = [-a2, -a1, a1, a2]
            self.w = [w2, w1, w1, w2]
        self.surfaces = [Surface(), Surface(), Surface(), Surface()]
        self.init_surfaces()
        self.pc2d = []
        for i in range(npc):
            for j in range(npc):
                ksi = self.pc[j]
                eta = self.pc[i]
                self.pc2d.append((ksi, eta))

    def pochodneKsi(self):
        for i in range(4):
            for p, (ksi, eta) in enumerate(self.pc2d):
                self.dN_dKsi[i][p] = dNx_po_dksi(eta, i)
        return self.dN_dKsi

    def pochodneEta(self):
        for i in range(4):
            for p, (ksi, eta) in enumerate(self.pc2d):
                self.dN_dEta[i][p] = dNx_po_deta(ksi, i)
        return self.dN_dEta
    def init_surfaces(self):
        g = Gauss(1)

        # BOK 0 → eta = -1
        for ksi in g.x:
            eta = -1
            self.surface_punkt(0, ksi, eta)

        # BOK 1 → ksi = 1
        for eta in g.x:
            ksi = 1
            self.surface_punkt(1, ksi, eta)

        # BOK 2 → eta = 1
        for ksi in g.x:
            eta = 1
            self.surface_punkt(2, ksi, eta)

        # BOK 3 → ksi = -1
        for eta in g.x:
            ksi = -1
            self.surface_punkt(3, ksi, eta)


    def surface_punkt(self, index, ksi, eta):
        N1 = 0.25 * (1 - ksi) * (1 - eta)
        N2 = 0.25 * (1 + ksi) * (1 - eta)
        N3 = 0.25 * (1 + ksi) * (1 + eta)
        N4 = 0.25 * (1 - ksi) * (1 + eta)
        self.surfaces[index].pc.append((ksi, eta))
        self.surfaces[index].N.append([N1, N2, N3, N4])


class Jakobian:
    def __init__(self, npc, element, grid, el_univ, data):
        self.npc = npc
        self.npc2d = npc * npc
        self.element = element
        self.grid = grid
        self.el_univ = el_univ
        self.data = data
        self.J = []
        self.detJ = []
        self.dNpodx = []
        self.dNpody = []
        self.H = []

    def obliczanieWartosci(self):
        dxdksi, dxdeta, dydksi, dydeta = [], [], [], []

        for j in range(self.npc2d):
            dx_dksi = dx_deta = dy_dksi = dy_deta = 0.0
            for i in range(4):
                node_id = self.element.nodeArray[i] - 1
                x = self.grid.nodesArray[node_id].x
                y = self.grid.nodesArray[node_id].y
                dx_dksi += self.el_univ.dN_dKsi[i][j] * x
                dx_deta += self.el_univ.dN_dEta[i][j] * x
                dy_dksi += self.el_univ.dN_dKsi[i][j] * y
                dy_deta += self.el_univ.dN_dEta[i][j] * y
            dxdksi.append(dx_dksi)
            dxdeta.append(dx_deta)
            dydksi.append(dy_dksi)
            dydeta.append(dy_deta)

        for i in range(self.npc2d):
            J = [[dxdksi[i], dxdeta[i]], [dydksi[i], dydeta[i]]]
            self.J.append(J)
            det = J[0][0] * J[1][1] - J[0][1] * J[1][0]
            self.detJ.append(det)

        for j in range(self.npc2d):
            odw_det = 1 / self.detJ[j]
            dNpodx1pc, dNpody1pc = [], []
            for i in range(4):
                dNdx = odw_det * (self.el_univ.dN_dKsi[i][j] * dydeta[j] - self.el_univ.dN_dEta[i][j] * dydksi[j])
                dNdy = odw_det * (-self.el_univ.dN_dKsi[i][j] * dxdeta[j] + self.el_univ.dN_dEta[i][j] * dxdksi[j])
                dNpodx1pc.append(dNdx)
                dNpody1pc.append(dNdy)
            self.dNpodx.append(dNpodx1pc)
            self.dNpody.append(dNpody1pc)

        H = [[0.0 for _ in range(4)] for _ in range(4)]
        for j in range(self.npc2d):
            detJ_local = self.detJ[j]
            w_ksi = self.el_univ.w[j % self.npc]
            w_eta = self.el_univ.w[j // self.npc]
            waga = w_ksi * w_eta
            for a in range(4):
                for b in range(4):
                    termX = self.dNpodx[j][a] * self.dNpodx[j][b]
                    termY = self.dNpody[j][a] * self.dNpody[j][b]
                    H[a][b] += self.data.Conductivity * (termX + termY) * detJ_local * waga
        self.H = H
        self.element.H = H


class UkladRownan:
    def __init__(self, grid1):
        self.grid = grid1
        self.hGlobalna = []

    def agregacja2d(self):
        H = [[0.0 for _ in range(self.grid.numNode)] for _ in range(self.grid.numNode)]
        for element in self.grid.elementsArray:
            hLokalna = element.H
            id = [n - 1 for n in element.nodeArray]
            for j in range(4):
                for k in range(4):
                    H[id[j]][id[k]] += hLokalna[j][k]
        self.hGlobalna = H

class HbcGenerator:
    def __init__(self, element, grid, data, el_univ):
        self.element = element
        self.grid = grid
        self.data = data
        self.el_univ = el_univ
        self.Hbc = [[0.0 for _ in range(4)] for _ in range(4)]
        self.w_1d = Gauss(1).w

    def dlugosc_boku(self, n1, n2):
        x1 = self.grid.nodesArray[n1 - 1].x
        y1 = self.grid.nodesArray[n1 - 1].y
        x2 = self.grid.nodesArray[n2 - 1].x
        y2 = self.grid.nodesArray[n2 - 1].y
        return math.sqrt((x2 - x1)**2 + (y2 - y1)**2)

    def compute(self):
        nodes = self.element.nodeArray

        boki = [
            (nodes[0], nodes[1], 0),
            (nodes[1], nodes[2], 1),
            (nodes[2], nodes[3], 2),
            (nodes[3], nodes[0], 3)
        ]

        for (n1, n2, s_id) in boki:
            if not (self.grid.nodesArray[n1-1].BC and self.grid.nodesArray[n2-1].BC):
                continue

            L = self.dlugosc_boku(n1, n2)
            detJ = L / 2.0

            surface = self.el_univ.surfaces[s_id]

            for p, N in enumerate(surface.N):
                w = self.w_1d[p]
                for i in range(4):
                    for j in range(4):
                        self.Hbc[i][j] += self.data.Alfa * N[i] * N[j] * w * detJ

        self.element.Hbc = self.Hbc
        return self.Hbc

