#! /usr/bin/python

# 6ta Practica Laboratorio 
# Complementos Matematicos I
# Ejemplo parseo argumentos

import argparse
import matplotlib.pyplot as plt
import numpy as np
from random import random, uniform
from math import sqrt


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

  def __init__(self, grafo, iters, refresh, temp, c1, c2, cooling, grav, verbose=False):
    """
    Parámetros:
    grafo: grafo en formato lista
    iters: cantidad de iteraciones a realizar
    refresh: cada cuántas iteraciones graficar. Si su valor es cero, entonces debe graficarse solo al final.
    temp: temperatura del sistema
    c1: constante de repulsión
    c2: constante de atracción
    cooling: Constante de enfriamiento
    verbose: si está encendido, activa los comentarios
    """

    # Guardo el grafo
    self.grafo = grafo
    self.vertices = grafo[0]
    self.aristas = grafo[1]

    # Inicializo estado
    self.posicion_x = {}
    self.posicion_y = {}
    self.accum_x = {}
    self.accum_y = {}
    self.largo = 1000 #Tamaño del frame

    # Guardo opciones
    self.iters = iters
    self.verbose = verbose
    self.temp = temp
    self.refresh = refresh
    self.c1 = c1
    self.c2 = c2
    #Guardamos las constantes K para no realizar demasiadas raices cuadradas en el algoritmo
    self.k1 = self.c1 * sqrt((self.largo)**2 / len(self.vertices))
    self.k2 = self.c2 * sqrt((self.largo)**2 / len(self.vertices))
    self.epsilon = 1 #Minima distancia entre vertices
    self.cooling = cooling
    self.gravedad = grav

  def inicializar_posiciones (self):
    for vertice in self.vertices:
      self.posicion_x[vertice] = uniform (0, self.largo)
      self.posicion_y[vertice] = uniform (0, self.largo)
  
  def inicializar_accum (self):
    for vertice in self.vertices:
      self.accum_x[vertice] = 0
      self.accum_y[vertice] = 0

  def actualizar_posicion (self):
    for vertice in self.vertices:
      fx = self.accum_x[vertice]
      fy = self.accum_y[vertice]
      modF = sqrt (fx**2 + fy**2)
      if modF > self.temp:
        fx = (fx / modF) * self.temp
        fy = (fy / modF) * self.temp
        self.accum_x[vertice] = fx
        self.accum_y[vertice] = fy

      nuevaPosicionX = self.posicion_x[vertice] + self.accum_x[vertice]
      nuevaPosicionY = self.posicion_y[vertice] + self.accum_y[vertice]
      
      if nuevaPosicionX > self.largo:
        self.posicion_x[vertice] = self.largo
      elif nuevaPosicionX < 0:
        self.posicion_x[vertice] = 0
      else: self.posicion_x[vertice] = nuevaPosicionX
      
      if nuevaPosicionY > self.largo:
        self.posicion_y[vertice] = self.largo
      elif nuevaPosicionY < 0:
        self.posicion_y[vertice] = 0
      else: self.posicion_y[vertice] = nuevaPosicionY
  
  def actualizar_temperatura (self):
    self.temp = self.cooling * self.temp
  
  def distancia (self, v1, v2):
    dist = sqrt((self.posicion_x[v2] - self.posicion_x[v1])**2 + (self.posicion_y[v2] - self.posicion_y[v1])**2)
    return dist
  
  def f_repulsion (self, dist):
    f = self.k1**2 / dist
    return f
  
  def f_atraccion (self, dist):
    f = dist**2 / self.k2
    return f

  def divison_por_cero (self, dist, v1, v2):
    while (dist < self.epsilon):
      fRandom = random()
      self.posicion_x[v1] += fRandom
      self.posicion_y[v1] += fRandom
      self.posicion_x[v2] -= fRandom
      self.posicion_y[v2] -= fRandom
      dist = self.distancia(v1, v2)
    return dist

  def calcular_gravedad (self):
    centro = self.largo/2
    for v in self.vertices:
      dist = sqrt ((self.posicion_x[v] - centro)**2 + (self.posicion_y[v] - centro)**2)
      #Caso division por cero
      while (dist < self.epsilon):
        fRandom = random()
        self.posicion_x[v] += fRandom
        self.posicion_y[v] += fRandom
        dist = sqrt ((self.posicion_x[v] - centro)**2 + (self.posicion_y[v] - centro)**2)

      fx = ((self.gravedad * (self.posicion_x[v] - centro)) / dist)
      fy = ((self.gravedad * (self.posicion_y[v] - centro)) / dist)
      self.accum_x[v] -= fx
      self.accum_y[v] -= fy
  
  def calcular_f_atrac (self):
    for [v1, v2] in self.aristas:
      distancia = self.distancia(v1, v2)
      # Caso de division por cero
      self.divison_por_cero(distancia, v1, v2)

      modF = self.f_atraccion (distancia)
      fx = (modF * (self.posicion_x[v2]-self.posicion_x[v1])) / distancia
      fy = (modF * (self.posicion_y[v2]-self.posicion_y[v1])) / distancia
      self.accum_x[v1] += fx
      self.accum_y[v1] += fy
      self.accum_x[v2] -= fx
      self.accum_y[v2] -= fy

  def calcular_f_rep (self):
    for v1 in self.vertices:
      for v2 in self.vertices:
        if v1 != v2:
          distancia = self.distancia(v1, v2)
          # Caso de divison por cero
          self.divison_por_cero(distancia, v1, v2)
          
          modF = self.f_repulsion (distancia)
          fx = (modF * (self.posicion_x[v2]-self.posicion_x[v1])) / distancia
          fy = (modF * (self.posicion_y[v2]-self.posicion_y[v1])) / distancia
          self.accum_x[v1] -= fx
          self.accum_y[v1] -= fy
          self.accum_x[v2] += fx
          self.accum_y[v2] += fy
  
  def step(self):
    self.inicializar_accum()
    self.calcular_f_atrac()
    self.calcular_f_rep()
    self.calcular_gravedad()
    self.actualizar_posicion()
    self.actualizar_temperatura()
  
  def dibujar (self):
    plt.pause(0.005)
    plt.clf()
    axes = plt.gca()
    axes.set_xlim([0, self.largo])
    axes.set_ylim([0, self.largo])
    plt.scatter (self.posicion_x.values(), self.posicion_y.values())
    for arista in self.aristas:
      vertice1 = arista[0]
      vertice2 = arista[1]
      plt.plot((self.posicion_x[vertice1], self.posicion_x[vertice2]), (self.posicion_y[vertice1], self.posicion_y[vertice2]))

  
  def layout(self):
    """
    Aplica el algoritmo de Fruchtermann-Reingold para obtener (y mostrar)
    un layout
    """
    self.inicializar_posiciones()
    plt.ion()
    for iter in range(self.iters):
      self.step()
      if self.refresh != 0 and (iter % self.refresh == 0):
        self.dibujar()
      elif self.refresh == 0:
        self.dibujar()
    plt.ioff()
    plt.show()

def main():
  # Definimos los argumentos de linea de comando que aceptamos
  parser = argparse.ArgumentParser()

  # Verbosidad, opcional, False por defecto
  parser.add_argument(
    '-v', '--verbose',
    action='store_true',
    help='Muestra mas informacion al correr el programa'
  )
  # Cantidad de iteraciones, opcional, 50 por defecto
  parser.add_argument(
    '--iters',
    type=int,
    help='Cantidad de iteraciones a efectuar',
    default=10
  )
  # Cantidad de refreshes, opcional, 1 por defecto
  parser.add_argument(
    '--refresh',
    type=int,
    help='Cada cuantas iteraciones graficar',
    default=1
  )
  # Temperatura inicial
  parser.add_argument(
    '--temp',
    type=float,
    help='Temperatura inicial',
    default=100.0
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
    default = 2
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
    default = 3
  )
  # Archivo del cual leer el grafo
  parser.add_argument(
    'file_name',
    help='Archivo del cual leer el grafo a dibujar'
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
    args.verbose
  )

  # Ejecutamos el layout
  layout_gr.layout()
  return


if __name__ == '__main__':
    main()