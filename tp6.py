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

  def __init__(self, grafo, iters, refresh, c1, c2, verbose=False):
    """
    Parámetros:
    grafo: grafo en formato lista
    iters: cantidad de iteraciones a realizar
    refresh: cada cuántas iteraciones graficar. Si su valor es cero, entonces debe graficarse solo al final.
    c1: constante de repulsión
    c2: constante de atracción
    verbose: si está encendido, activa los comentarios
    """

    # Guardo el grafo
    self.grafo = grafo
    self.vertices = grafo[0]
    self.aristas = grafo[1]

    # Inicializo estado
    # Completar
    self.posicion_x = {}
    self.posicion_y = {}
    self.accum_x = {}
    self.accum_y = {}
    self.largo = 100
    self.fuerzas_atractivas = {}
    self.fuerzas_repulsivas = {}

    # Guardo opciones
    self.iters = iters
    self.verbose = verbose
    # TODO: faltan opciones TEMPERATURA
    self.refresh = refresh
    self.c1 = c1
    self.c2 = c2

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
      #CUIDADO CON LOS BORDES DE LA VENTANA
      self.posicion_x[vertice] += self.accum_x[vertice]
      self.posicion_y[vertice] += self.accum_y[vertice]
  
  def distancia (self, v1, v2):
    dist = sqrt((self.posicion_x[v2] - self.posicion_x[v1])**2 + (self.posicion_y[v2] - self.posicion_y[v1])**2)
    return dist
  
  def f_repulsion (self, dist):
    k = self.c1 * sqrt((self.largo)**2 / len(self.vertices))
    f = k**2 / dist
    return f
  
  def f_atraccion (self, dist):
    k = self.c2 * sqrt ((self.largo)**2 / len(self.vertices))
    f = dist**2 / k
    return f

  def calcular_f_atrac (self):
    for [v1, v2] in self.aristas:
      distancia = self.distancia(v1, v2)
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
    self.actualizar_posicion()
  
  def dibujar (self):
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
    self.dibujar()
    plt.show()
    for iter in range(self.iters):
      self.step()
      self.dibujar()
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
  # Temperatura inicial
  parser.add_argument(
    '--temp',
    type=float,
    help='Temperatura inicial',
    default=100.0
  )
  # Archivo del cual leer el grafo
  parser.add_argument(
    'file_name',
    help='Archivo del cual leer el grafo a dibujar'
  )

  args = parser.parse_args()

  # Descomentar abajo para ver funcionamiento de argparse
  # print args.verbose
  # print args.iters    
  # print args.file_name
  # print args.temp
  # return

  # Creamos nuestro objeto LayoutGraph
  layout_gr = LayoutGraph(
    grafo = leer_archivo (args.file_name),
    iters = args.iters,
    refresh = 1,
    c1=0.1,
    c2=5.0,
    verbose=False
  )

  # Ejecutamos el layout
  layout_gr.layout()
  return


if __name__ == '__main__':
    main()