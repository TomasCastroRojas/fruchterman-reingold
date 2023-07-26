#! /usr/bin/python

# Trabajo Practico Final
# Complementos Matematicos I
# Blas Barbagelata - Tomas Castro Rojas

# Linea de comando para ejecutar el algoritmo:
# python3 Barbagelata_CastroRojas.py file_name 
# Banderas opcionales (sin corchetes):
# [-h] [-v] [--iters ITERS] [--refresh REFRESH] [--temp TEMP] [--c1 C1] [--c2 C2] 
# [--cooling COOLING] [--grav GRAV] [--largo LARGO]
# Si el argumento refresh es 0, imprime solo el estado final del grafo.

import argparse
import matplotlib.pyplot as plt
from random import random, uniform
from math import sqrt
import numpy as np

Vertice : type = str
Vector2D : type = tuple[float, float]
Posicion : type = dict[Vertice, Vector2D]
Fuerza : type = dict[Vertice, Vector2D]
Grafo : type = tuple[list[str], list[str]]

# Lee un archivo que contiene un grafo, devuelve el grafo en forma G = (V,E)
def leer_archivo (filename: str) -> Grafo: 
    with open(filename, 'r') as file:
        lineas = file.readlines()
    n = int(lineas[0])
    vertices = [v.strip() for v in lineas[1:n+1]]
    aristas = [e.split() for e in lineas[n+1:]]
    return vertices, aristas

class LayoutGraph:

    def __init__(self, grafo: tuple[list[str],list[str]], iters: int, refresh: int,
                       temp: float, c1: float, c2: float, cooling: float, grav: float, largo: int, verbose=False):
        """
        Parámetros:
        grafo: grafo en formato lista
        iters: cantidad de iteraciones a realizar
        refresh: cada cuántas iteraciones graficar. Si su valor es cero, entonces debe graficarse solo al final.
        temp: temperatura del sistema
        c1: constante de repulsión
        c2: constante de atracción
        cooling: Constante de enfriamiento
        grav: constante de gravedad
        largo: Tamaño del frame
        verbose: si está encendido, activa los comentarios
        """

        # Guardo el grafo
        self.vertices = grafo[0]
        self.aristas = grafo[1]

        # Inicializo estado
        self.posiciones : Posicion = {v: np.random.rand(2)*largo for v in self.vertices}
        self.acum : Fuerza = {}
        self.largo = largo #Tamaño del frame

        # Guardo opciones
        self.iters = iters
        self.verbose = verbose
        self.temp = temp
        self.refresh = refresh
        self.c1 = c1
        self.c2 = c2
        self.cooling = cooling
        self.gravedad = grav

        # Guardamos las constantes K para no realizar demasiadas raices cuadradas en el algoritmo
        ratio = sqrt((self.largo)**2 / len(self.vertices))
        self.k1 : float = self.c1 * ratio
        self.k2 : float = self.c2 * ratio
        self.epsilon : float = 1 #Minima distancia entre vertices
        self.center : Vector2D = np.array((self.largo/2, self.largo/2))
    
    # Inicializa los acumuladores de fuerza en 0
    def inicializar_acum (self):
        self.imprimir_mensaje(f"--Inicializando acumuladores de fuerza--")
        self.acum = {v: np.zeros(2) for v in self.vertices}

    # Si la distancia entre dos vertices es menor a la minima, aplica fuerzas
    # aleatorias a los vertices hasta que supera la distancia minima
    def divison_por_cero (self, dist: float, v1: Vertice, v2: Vertice):
        e = self.posiciones[v2] - self.posiciones[v1]
        while (dist < self.epsilon):
            self.imprimir_mensaje(f"--Distancia menor a la minima. Aplicando fuerza aleatoria--")

            fRandom = np.random.rand(2)
            self.posiciones[v1] += fRandom
            self.posiciones[v2] -= fRandom
            e = self.posiciones[v2] - self.posiciones[v1]
            dist = np.linalg.norm(e)
        return dist, e

    # Calcula las fuerzas de atraccion entre las aristas
    def calcular_f_atrac (self):
        self.imprimir_mensaje(f"--Calculando fuerzas de atraccion--")

        for v1, v2 in self.aristas:
            e = self.posiciones[v2] - self.posiciones[v1]
            dist = np.linalg.norm(e)
            dist, e = self.divison_por_cero(dist, v1, v2)

            moduloF = dist / self.k2
            f = moduloF * e
            self.acum[v1] += f
            self.acum[v2] -= f

        self.imprimir_mensaje(f"--Fin calculo fuerzas de atraccion--")

    # Calcula las fuerzas de repulsion entre vertices
    def calcular_f_rep (self):
        self.imprimir_mensaje(f"--Calculando fuerzas de repulsion--")

        for v1 in self.vertices:
            for v2 in self.vertices:
                if v1 != v2:
                    e = self.posiciones[v2] - self.posiciones[v1]
                    dist = np.linalg.norm(e)
                    dist, e = self.divison_por_cero(dist, v1, v2)
                    
                    moduloF = (self.k1 / dist)**2
                    f = moduloF * e
                    self.acum[v1] -= f
                    self.acum[v2] += f

        self.imprimir_mensaje(f"--Fin calculo fuerzas de repulsion--")
    
    # Calcula la fuerza de gravedad que se ejerce sobre cada vertice
    def calcular_gravedad (self):
        self.imprimir_mensaje(f"--Calculando fuerzas de gravedad--")

        for v in self.vertices:
            e = self.posiciones[v] - self.center
            dist = np.linalg.norm(e)
            #Caso division por cero
            while (dist < self.epsilon):
                fRandom = np.random.rand(2)
                self.posiciones[v] += fRandom
                e = self.posiciones[v] - self.center
                dist = np.linalg.norm(e)

            f = (self.gravedad/dist) * e
            self.acum[v] -= f

        self.imprimir_mensaje(f"--Fin calculo fuerzas de gravedad--")

    # Actualiza las posiciones de los vertices
    def actualizar_posicion (self):
        self.imprimir_mensaje(f"--Actualizando posiciones de los vertices--")

        for v in self.vertices:
            f = self.acum[v]
            moduloF = np.linalg.norm(f)
            if moduloF > self.temp:
                self.acum[v] *= self.temp / moduloF

            self.posiciones[v] += self.acum[v]
            
            np.clip(self.posiciones[v], 0, self.largo, out=self.posiciones[v])
          
        self.imprimir_mensaje(f"--Nuevas posiciones calculadas--")
    
    # Actualiza la temperatura
    def actualizar_temperatura (self):
        self.imprimir_mensaje(f"--Actualizando temperatura--")
        self.temp = self.cooling * self.temp

    # Rutinas a aplicar en cada iteracion del algoritmo
    def step(self):
        self.inicializar_acum()
        self.calcular_f_atrac()
        self.calcular_f_rep()
        self.calcular_gravedad()
        self.actualizar_posicion()
        self.actualizar_temperatura()
    
    # Dibuja los vertices y aristas
    def actualizar_plot (self):
        plt.clf()
        plt.xlim(0, self.largo)
        plt.ylim(0, self.largo)

        for v1, v2 in self.aristas:
            p1 = self.posiciones[v1]
            p2 = self.posiciones[v2]
            plt.plot([p1[0], p2[0]], [p1[1], p2[1]])
        
        xcoord = [p[0] for p in self.posiciones.values()]
        ycoord = [p[1] for p in self.posiciones.values()]

        plt.scatter(xcoord, ycoord)

        plt.show()
        plt.pause(0.005)

    # Si la opcion verbose esta activada, imprime un mensaje
    def imprimir_mensaje(self, mensaje: str):
        if self.verbose:
            print(mensaje)

    def layout(self):
        """
        Aplica el algoritmo de Fruchtermann-Reingold para obtener (y mostrar)
        un layout
        """
        plt.ion()
        self.actualizar_plot()
        for iter in range(self.iters):
            self.imprimir_mensaje(f"\n--Iteracion {iter}--")

            self.step()
            if self.refresh != 0 and (iter % self.refresh == 0):
                self.imprimir_mensaje(f"--Dibuja el grafo en iteracion {iter}--")

                self.actualizar_plot()
        self.imprimir_mensaje(f"--Fin del algoritmo--")
        plt.ioff()
        self.actualizar_plot()

def main():
    # Definimos los argumentos de linea de comando que aceptamos
    parser = argparse.ArgumentParser()

    # Verbosidad
    parser.add_argument(
        '-v', '--verbose',
        action = 'store_true',
        help = 'Muestra mas informacion al correr el programa'
    )
    # Cantidad de iteraciones
    parser.add_argument(
        '--iters',
        type = int,
        help = 'Cantidad de iteraciones a efectuar',
        default = 100
    )
    # Cantidad de refreshes 0 por defecto
    parser.add_argument(
        '--refresh',
        type = int,
        help = 'Cada cuantas iteraciones graficar',
        default = 5
    )
    # Temperatura inicial
    parser.add_argument(
        '--temp',
        type = float,
        help = 'Temperatura inicial',
        default = 100.0
    )
    # Constante de repulsion
    parser.add_argument(
        '--c1',
        type = float,
        help = 'Constante de repulsion',
        default = 0.1
    )
    # Constante de atraccion
    parser.add_argument(
        '--c2',
        type = float,
        help = 'Constante de atraccion',
        default = 3
    )
    # Constante de enfriamiento
    parser.add_argument(
        '--cooling',
        type = float,
        help = 'Constante de enfriamiento',
        default = 0.95
    )
    # Constante de gravedad
    parser.add_argument(
        '--grav',
        type = float,
        help = 'Constante de gravedad',
        default = 2.5
    )
    # Largo del frame
    parser.add_argument(
        '--largo',
        type = int,
        help = 'Tamaño del frame',
        default = 500
    )
    # Archivo del cual leer el grafo
    parser.add_argument(
        'file_name',
        help = 'Archivo del cual leer el grafo a dibujar'
    )

    args = parser.parse_args()
    graph = leer_archivo(args.file_name)

    # Creamos nuestro objeto LayoutGraph
    layout_gr = LayoutGraph(
        graph,
        args.iters,
        args.refresh,
        args.temp,
        args.c1,
        args.c2,
        args.cooling,
        args.grav,
        args.largo,
        args.verbose
    )

    # Ejecutamos el layout
    layout_gr.layout()
    return


if __name__ == '__main__':
    main()