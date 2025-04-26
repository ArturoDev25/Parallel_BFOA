# Parallel_BFOA

Este proyecto implementa y mejora el algoritmo **Bacterial Foraging Optimization Algorithm (BFOA)** para el **alineamiento de secuencias biológicas** usando **matrices BLOSUM**, procesamiento paralelo y técnicas de optimización bioinspiradas.

## 📌 Objetivo

Alinear múltiples secuencias de proteínas de forma óptima, maximizando el puntaje de similitud (BLOSUM) y la interacción entre soluciones, utilizando BFOA como algoritmo base y evaluando mejoras aplicadas al mismo.

---

## 📁 Estructura del Proyecto

- `bacteria.py`: Versión original del algoritmo BFOA.
- `bacteria_Mejorada.py`: Versión mejorada con tumbo inteligente y penalización por gaps.
- `parallel_BFOA.py`: Script para ejecutar 30 corridas del algoritmo original.
- `parallel_BFOA_Mejorada.py`: Script para ejecutar 30 corridas del algoritmo mejorado.
- `Comparar.ipynb`: Análisis comparativo con visualizaciones entre ambos algoritmos.
- `evaluadorBlosum.py`: Calcula puntajes de similitud usando la matriz BLOSUM.
- `fastaReader.py`: Lector de archivos multiFASTA.
- `Data.zip`: Secuencias de entrada en formato FASTA.
- `Resultados_BFOA.csv`: Resultados del algoritmo original.
- `Resultados_BFOA_Mejorado.csv`: Resultados del algoritmo mejorado.

---

## 🚀 Mejoras implementadas

### 1. Tumbo Inteligente
En lugar de insertar gaps aleatoriamente, la mejora detecta regiones de alta variabilidad en las secuencias y aplica los gaps donde generan mayor beneficio biológico.

### 2. Penalización por Gaps
Cada gap insertado resta puntuación al fitness, evitando soluciones artificiales que maximizan el score BLOSUM a costa de introducir demasiados espacios.

---

## 📊 Resultados

- Las soluciones mejoradas son más **estables** y **biológicamente razonables**.
- Se reduce el uso excesivo de gaps sin sacrificar significativamente el fitness global.
- Se mantiene una interacción efectiva entre bacterias.

> Consulta el archivo `Comparar.ipynb` para visualizar histogramas, evolución métrica y análisis de correlación.

---

## 🧪 Requisitos

- Python 3.8+
- pandas
- matplotlib
- seaborn
- numpy
- scikit-learn
