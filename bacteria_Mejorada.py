from multiprocessing import Manager, Pool
from evaluadorBlosum import evaluadorBlosum
import numpy as np
import random
import copy

class bacteria_mejorada:

    def __init__(self, numBacterias):
        manager = Manager()
        self.blosumScore = manager.list(range(numBacterias))
        self.tablaAtract = manager.list(range(numBacterias))
        self.tablaRepel = manager.list(range(numBacterias))
        self.tablaInteraction = manager.list(range(numBacterias))
        self.tablaFitness = manager.list(range(numBacterias))
        self.granListaPares = manager.list(range(numBacterias))
        self.NFE = manager.list(range(numBacterias))
        self.contadorGaps = manager.list([0] * numBacterias)

    def resetListas(self, numBacterias):
        manager = Manager()
        self.__init__(numBacterias)

    def cuadra(self, numSec, poblacion):
        for i in range(len(poblacion)):
            bacterTmp = list(poblacion[i])[:numSec]
            maxLen = max(len(seq) for seq in bacterTmp)
            for t in range(numSec):
                gap_count = maxLen - len(bacterTmp[t])
                if gap_count > 0:
                    bacterTmp[t].extend(["-"] * gap_count)
            poblacion[i] = tuple(bacterTmp)

    def tumbo_inteligente(self, numSec, poblacion, numGaps=200):
        for i in range(len(poblacion)):
            bacterTmp = list(poblacion[i])
            variability = self.calculaVariabilidad(bacterTmp)
            for _ in range(numGaps):
                seq_idx = random.randint(0, numSec - 1)
                pos = self.escogePosicionConAltaVariabilidad(variability)
                bacterTmp[seq_idx].insert(pos, "-")
                self.contadorGaps[i] += 1
            poblacion[i] = tuple(bacterTmp)

    def calculaVariabilidad(self, secuencias):
        numCols = max(len(seq) for seq in secuencias)
        variabilidad = []
        for col in range(numCols):
            columna = [seq[col] if col < len(seq) else "-" for seq in secuencias]
            freqs = {c: columna.count(c) for c in set(columna)}
            probs = [f / len(columna) for f in freqs.values()]
            entropia = -sum(p * np.log2(p) for p in probs if p > 0)
            variabilidad.append(entropia)
        return variabilidad

    def escogePosicionConAltaVariabilidad(self, variabilidad):
        total = sum(variabilidad)
        if total == 0:
            return random.randint(0, len(variabilidad) - 1)
        probs = [v / total for v in variabilidad]
        return np.random.choice(len(variabilidad), p=probs)

    def penalizaPorGaps(self, idx):
        return self.contadorGaps[idx] * 20

    def creaTablaFitness(self):
        for i in range(len(self.tablaInteraction)):
            valorBlsm = self.blosumScore[i]
            valorInteract = self.tablaInteraction[i]
            penalizacion = self.penalizaPorGaps(i)
            self.tablaFitness[i] = valorBlsm + valorInteract - penalizacion

    def creaGranListaPares(self, poblacion):
        for i in range(len(poblacion)):
            pares = []
            bacterTmp = list(poblacion[i])
            for j in range(len(bacterTmp[0])):
                columna = [seq[j] for seq in bacterTmp]
                pares.extend(self.obtener_pares_unicos(columna))
            self.granListaPares[i] = pares

    def obtener_pares_unicos(self, columna):
        pares = set()
        for i in range(len(columna)):
            for j in range(i+1, len(columna)):
                par = tuple(sorted([columna[i], columna[j]]))
                pares.add(par)
        return list(pares)

    def evaluaFila(self, fila, idx):
        score = sum(evaluadorBlosum().getScore(a, b) for a, b in fila if a != '-' and b != '-')
        self.blosumScore[idx] = score

    def evaluaBlosum(self):
        with Pool() as pool:
            args = [(copy.deepcopy(self.granListaPares[i]), i) for i in range(len(self.granListaPares))]
            pool.starmap(self.evaluaFila, args)

    def compute_diff(self, args):
        indexBacteria, otherScore, selfScores, d, w = args
        diff = min((selfScores[indexBacteria] - otherScore) ** 2, 500)
        self.NFE[indexBacteria] += 1
        return d * np.exp(w * diff)

    def compute_cell_interaction(self, indexBacteria, d, w, atracTrue):
        with Pool() as pool:
            args = [(indexBacteria, other, self.blosumScore, d, w) for other in self.blosumScore]
            results = pool.map(self.compute_diff, args)
        total = sum(results)
        if atracTrue:
            self.tablaAtract[indexBacteria] = total
        else:
            self.tablaRepel[indexBacteria] = total

    def creaTablasAtractRepel(self, poblacion, dAttr, wAttr, dRepel, wRepel):
        for idx in range(len(poblacion)):
            self.compute_cell_interaction(idx, dAttr, wAttr, True)
            self.compute_cell_interaction(idx, dRepel, wRepel, False)

    def creaTablaInteraction(self):
        for i in range(len(self.tablaAtract)):
            self.tablaInteraction[i] = self.tablaAtract[i] + self.tablaRepel[i]

    def getNFE(self):
        return sum(self.NFE)

    def obtieneBest(self, globalNFE):
        bestIdx = max(range(len(self.tablaFitness)), key=lambda i: self.tablaFitness[i])
        print(f"BEST #{bestIdx}: Fitness={self.tablaFitness[bestIdx]}, BLOSUM={self.blosumScore[bestIdx]}, Interaction={self.tablaInteraction[bestIdx]}, Gaps={self.contadorGaps[bestIdx]}, NFE={globalNFE}")
        return bestIdx, self.tablaFitness[bestIdx]

    def replaceWorst(self, poblacion, bestIdx):
        worstIdx = min(range(len(self.tablaFitness)), key=lambda i: self.tablaFitness[i])
        poblacion[worstIdx] = copy.deepcopy(poblacion[bestIdx])
        self.contadorGaps[worstIdx] = self.contadorGaps[bestIdx]
