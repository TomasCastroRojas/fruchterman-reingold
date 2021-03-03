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


# Lee un archivo que contiene un grafo, devuelve el grafo en forma G = [V,E]
def leer_archivo (nombreArchivo):
  grafoArchivo = open (nombreArchivo, "r")
  grafo = [[], []]
  cantVertices = int(grafoArchivo.readline().rstrip("\n"))
  lineas = grafoArchivo.readlines()
  grafoArchivo.close()
  cantLineas = len(lineas)
  #Agrego los vertices
  for vertices in range(cantVertices):
    grafo[0].append(lineas[vertices].rstrip("\n"))
  #Agrego las aristas
  for aristas in range(cantVertices,cantLineas):
    grafo[1].append(lineas[aristas].split())
  return grafo

class LayoutGraph:

  def __init__(self, grafo, iters, refresh, temp, c1, c2, cooling, grav, largo, verbose=False):
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
    self.grafo = grafo
    self.vertices = grafo[0]
    self.aristas = grafo[1]

    # Inicializo estado
    self.posicionX = {}
    self.posicionY = {}
    self.acumX = {}
    self.acumY = {}
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
    self.k1 = self.c1 * ratio
    self.k2 = self.c2 * ratio
    self.epsilon = 1 #Minima distancia entre vertices

  # Inicializa las posiciones de los vertices de forma aleatoria
  def inicializar_posiciones (self):
    self.imprimir_mensaje(f"--Inicializando posiciones de los vertices--")
    for vertice in self.vertices:
      self.posicionX[vertice] = uniform (0, self.largo)
      self.posicionY[vertice] = uniform (0, self.largo)
  
  # Inicializa los acumuladores de fuerza en 0
  def inicializar_acum (self):
    self.imprimir_mensaje(f"--Inicializando acumuladores de fuerza--")
    for vertice in self.vertices:
      self.acumX[vertice] = 0
      self.acumY[vertice] = 0
  
  # Calcula la distancia euclidiana entre dos vertices
  def distancia (self, v1, v2):
    dist = sqrt((self.posicionX[v2] - self.posicionX[v1])**2 + (self.posicionY[v2] - self.posicionY[v1])**2)
    return dist
  
  # Calcula la fuerza de atraccion
  def f_atraccion (self, dist):
    f = dist**2 / self.k2
    return f

  # Calcula la fuerza de repulsion
  def f_repulsion (self, dist):
    f = self.k1**2 / dist
    return f

  # Si la distancia entre dos vertices es menor a la minima, aplica fuerzas
  # aleatorias a los vertices hasta que supera la distancia minima
  def divison_por_cero (self, dist, v1, v2):
    while (dist < self.epsilon):
      self.imprimir_mensaje(f"--Distancia menor a la minima. Aplicando fuerza aleatoria--")
      fRandom = random()
      self.posicionX[v1] += fRandom
      self.posicionY[v1] += fRandom
      self.posicionX[v2] -= fRandom
      self.posicionY[v2] -= fRandom
      dist = self.distancia(v1, v2)
    return dist

  # Calcula las fuerzas de atraccion entre las aristas
  def calcular_f_atrac (self):
    self.imprimir_mensaje(f"--Calculando fuerzas de atraccion--")
    for [v1, v2] in self.aristas:
      dist = self.distancia(v1, v2)
      # Caso de division por cero
      dist = self.divison_por_cero(dist, v1, v2)

      moduloF = self.f_atraccion (dist)
      fx = (moduloF * (self.posicionX[v2]-self.posicionX[v1])) / dist
      fy = (moduloF * (self.posicionY[v2]-self.posicionY[v1])) / dist
      self.acumX[v1] += fx
      self.acumY[v1] += fy
      self.acumX[v2] -= fx
      self.acumY[v2] -= fy
    self.imprimir_mensaje(f"--Fin calculo fuerzas de atraccion--")

  # Calcula las fuerzas de repulsion entre vertices
  def calcular_f_rep (self):
    self.imprimir_mensaje(f"--Calculando fuerzas de repulsion--")
    for v1 in self.vertices:
      for v2 in self.vertices:
        if v1 != v2:
          dist = self.distancia(v1, v2)
          # Caso de divison por cero
          dist = self.divison_por_cero(dist, v1, v2)
          
          moduloF = self.f_repulsion (dist)
          fx = (moduloF * (self.posicionX[v2]-self.posicionX[v1])) / dist
          fy = (moduloF * (self.posicionY[v2]-self.posicionY[v1])) / dist
          self.acumX[v1] -= fx
          self.acumY[v1] -= fy
          self.acumX[v2] += fx
          self.acumY[v2] += fy
    self.imprimir_mensaje(f"--Fin calculo fuerzas de repulsion--")
  
  # Calcula la fuerza de gravedad que se ejerce sobre cada vertice
  def calcular_gravedad (self):
    self.imprimir_mensaje(f"--Calculando fuerzas de gravedad--")
    centro = self.largo/2
    for v in self.vertices:
      dist = sqrt ((self.posicionX[v] - centro)**2 + (self.posicionY[v] - centro)**2)
      #Caso division por cero
      while (dist < self.epsilon):
        fRandom = random()
        self.posicionX[v] += fRandom
        self.posicionY[v] += fRandom
        dist = sqrt ((self.posicionX[v] - centro)**2 + (self.posicionY[v] - centro)**2)

      fx = ((self.gravedad * (self.posicionX[v] - centro)) / dist)
      fy = ((self.gravedad * (self.posicionY[v] - centro)) / dist)
      self.acumX[v] -= fx
      self.acumY[v] -= fy
    self.imprimir_mensaje(f"--Fin calculo fuerzas de gravedad--")

  # Actualiza las posiciones de los vertices
  def actualizar_posicion (self):
    self.imprimir_mensaje(f"--Actualizando posiciones de los vertices--")
    for vertice in self.vertices:
      fx = self.acumX[vertice]
      fy = self.acumY[vertice]
      moduloF = sqrt (fx**2 + fy**2)
      if moduloF > self.temp:
        fx = (fx / moduloF) * self.temp
        fy = (fy / moduloF) * self.temp
        self.acumX[vertice] = fx
        self.acumY[vertice] = fy

      nuevaPosicionX = self.posicionX[vertice] + self.acumX[vertice]
      nuevaPosicionY = self.posicionY[vertice] + self.acumY[vertice]
      
      if nuevaPosicionX > self.largo:
        self.imprimir_mensaje(f"--Vertice fuera de rango. Acomodando--")
        self.posicionX[vertice] = self.largo
      elif nuevaPosicionX < 0:
        self.imprimir_mensaje(f"--Vertice fuera de rango. Acomodando--")
        self.posicionX[vertice] = 0
      else: self.posicionX[vertice] = nuevaPosicionX
      
      if nuevaPosicionY > self.largo:
        self.imprimir_mensaje(f"--Vertice fuera de rango. Acomodando--")
        self.posicionY[vertice] = self.largo
      elif nuevaPosicionY < 0:
        self.imprimir_mensaje(f"--Vertice fuera de rango. Acomodando--")
        self.posicionY[vertice] = 0
      else: self.posicionY[vertice] = nuevaPosicionY
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
    plt.pause(0.005)
    plt.clf()
    axes = plt.gca()
    axes.set_xlim([0, self.largo])
    axes.set_ylim([0, self.largo])
    plt.scatter (self.posicionX.values(), self.posicionY.values())
    for arista in self.aristas:
      v1 = arista[0]
      v2 = arista[1]
      plt.plot((self.posicionX[v1], self.posicionX[v2]), (self.posicionY[v1], self.posicionY[v2]))

  # Si la opcion verbose esta activada, imprime un mensaje
  def imprimir_mensaje(self, mensaje):
    if self.verbose:
      print(mensaje)

  def layout(self):
    """
    Aplica el algoritmo de Fruchtermann-Reingold para obtener (y mostrar)
    un layout
    """
    self.inicializar_posiciones()
    plt.ion()
    self.actualizar_plot()
    for iter in range(self.iters):
      self.imprimir_mensaje(f"\n--Iteracion {iter}--")
      self.step()
      if self.refresh != 0 and (iter % self.refresh == 0):
        self.imprimir_mensaje(f"--Dibuja el grafo en iteracion {iter}--")
        self.actualizar_plot()
    plt.ioff()
    self.imprimir_mensaje(f"--Fin del algoritmo--")
    self.actualizar_plot()
    plt.show() # Imprime el estado final

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
    default = 50
  )
  # Cantidad de refreshes 0 por defecto
  parser.add_argument(
    '--refresh',
    type = int,
    help = 'Cada cuantas iteraciones graficar',
    default = 1
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

  # Creamos nuestro objeto LayoutGraph
  layout_gr = LayoutGraph(
    leer_archivo (args.file_name),
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