#!pip install gurobipy

import random
import networkx as nx
import gurobipy as gp
from gurobipy import GRB
import numpy as np
from datetime import datetime
import math
import matplotlib.pyplot as plt
import time
from queue import PriorityQueue
import csv
import os
import sys
# sys.path.append("./miniforge3/lib/python3.10/site-packages")
# import latextable
# from texttable import Texttable

maxDemand = 100
minDemand = 0
numMatching = 1
degree= 2
fixedDemand = 1000
thersholdDemand = 100

"""## Creating Input

"""

class Instance:

  def read_edges_from_file(filename):
    edges = []
    curFileName = filename
    with open(curFileName, 'r') as file:
        for line in file:
            node1, node2, weight = line.split()
            edges.append((node1, node2, weight))
    return edges

  def getRealDem(filename,numNodes):
    edges = Instance.read_edges_from_file(filename)
    nodes = set()
    n = 0
    for edge in edges:
        n = max(n,int(edge[0]))
        n = max(n,int(edge[1]))
    n = n+1
    dem = [[0] * n for _ in range(n)]
    for u, v, weight in edges:
      xU = int(u)
      xV = int(v)
      if (xU< n and xV < n):
        dem[xU][xV] = int(weight)
    return dem


  def createRing(myNumNodes):
    return nx.cycle_graph(myNumNodes)

  def createGrid(myNumNodes):
    dimension = int(math.sqrt(myNumNodes))
    return nx.grid_2d_graph(dimension, dimension)

  def notAnEdge(self,u,v):
    return not self.graph.has_edge(u, v)

  def printGraph(nowGraph, drawPos):
    drawColors = [nowGraph[u][v]['color'] for u,v in nowGraph.edges]
    drawAlphas = [nowGraph[u][v]['alpha'] for u,v in nowGraph.edges]
    drawWidth = [nowGraph[u][v]['width'] for u,v in nowGraph.edges]
    nx.draw_networkx_edges(nowGraph,pos = drawPos, edge_color=drawColors, alpha=drawAlphas, width=drawWidth)

  def setDemand(self,inPutDemand):
    self.demand = inPutDemand

  def twoDtoOneD(x,y,n):
    return x+y*n

  def Torus2D(myNumNodes):
    n = round(myNumNodes ** (1 / 2))
    G = nx.Graph()
    for x in range(n):
      for y in range(n):
          for direction in [[1, 0], [0, 1]]:
              xx, yy = (x + direction[0]) % n, (y + direction[1]) % n
              G.add_edge(Instance.twoDtoOneD(x, y, n), Instance.twoDtoOneD(xx, yy, n))
    return G

  def threeDtoOneD(x,y,z,n):
    return x+ y*n + z*(n**2)

  def Torus3D(myNumNodes):
    n = round(myNumNodes ** (1 / 3))
    G = nx.Graph()
    for x in range(n):
      for y in range(n):
        for z in range(n):
          for dx, dy, dz in [(1, 0, 0), (0, 1, 0), (0, 0, 1)]:
              xx, yy, zz = (x + dx) % n, (y + dy) % n, (z + dz) % n
              G.add_edge(Instance.threeDtoOneD(x, y, z, n), Instance.threeDtoOneD(xx, yy, zz,n))

    return G

  def createInfrastructure(myNumNodes, infType):
    if(infType == "ring"):
      return Instance.createRing(myNumNodes)
    if (infType == "2DTorus"):
      return Instance.Torus2D(myNumNodes)
    if (infType == "3DTorus"):
      return Instance.Torus3D(myNumNodes)

  def findClosetPow2(n):
    return 2 ** round(np.log2(n))

  def transformToPowerTwo(inst):
    currN = inst.numNodes
    twoN = findClosetPow2(currN)
    inst.numNodes = twoN
    return inst

  def createSuitDemand(dmenad_edges, n):
    dem = [[0.0] * n for _ in range(n)]
    for u, v, weight in dmenad_edges:
      xU = int(u)
      xV = int(v)
      if (xU< n and xV < n):
        dem[xU][xV] = float(weight)
    return dem

  def getDemand(self):
    return self.demand

  def getGraph(self):
    return self.graph

  def __init__(self, infType, demandType, demand_edges, graphNumNodes):
    self.numNodes = graphNumNodes
    self.demand = Instance.createSuitDemand(demand_edges,self.numNodes)
    self.MaxValue = self.numNodes +2
    self.graph = Instance.createInfrastructure(self.numNodes, infType)
    self.nodes = self.graph.nodes
    self.edges = self.graph.edges
    # if(demandType == "ringDem"): #Creating ring demand
    #   self.demand, self.curDrawGraph = Instance.createRingDemand(genValue,myNumNodes)
    # if(demandType == "random"): #Creating random demand
    #   self.demand = Instance.createRandomDemand(myNumNodes,genValue)
    # if(".txt" in demandType):
    #   self.demand = Instance.getRealDem(demandType, myNumNodes) #TODO
    # if(demandType == "zipf"): #Creating random demand
    #   self.demand = Instance.createZipfDemand(myNumNodes,genValue)
    #self.drawPos = nx.circular_layout(self.graph)

"""# Main to Run

## Solving for different modes

## Main running
"""

def read_weighted_edges_from_folder(folder_path):
    all_edges = []

    for filename in os.listdir(folder_path):
        file_path = os.path.join(folder_path, filename)
        if os.path.isfile(file_path):
            edges = []
            with open(file_path, mode='r') as file:
                csv_reader = csv.reader(file)
                [extra]= next(csv_reader)  # Skip the first row
                numNodes = int(extra.split()[1])
                for row in csv_reader:
                    u, v, weight = row
                    edges.append((int(u)-1, int(v)-1, float(weight)))
            all_edges.append([str(filename),numNodes,edges])
    return all_edges

def saveCSV(graphRes):
  resRes = []
  resRes.append(["Name", "$|N|$", "$|E|$", "LP", "FG", "MoD", "SC", "SD",  "Base"])
  for i in range(len(graphRes)):
    x = [graphRes[i][0]]
    x.append(str(graphRes[i][1]))
    x.append(str(graphRes[i][2]))
    for y in graphRes[i][3]:
      x.append(str(y))
    resRes.append(x)
  #print(resRes)
  resTime = []

  # for i in range(len(graphRes)):
  #   resTime.append(graphRes[i][1])
  nameCSV = "SuiteResult" + datetime.now().strftime('%Y-%m-%d-%H:%M:%S') + ".csv"
  np.savetxt(nameCSV, resRes, fmt="%s", delimiter=',')

all_edges = read_weighted_edges_from_folder("/Users/pourdamghani/Desktop/SparseSuite")
graphRes = []
for graph_Data in all_edges:
  graph_edges = graph_Data[2]
  graphNumNodes = graph_Data[1]
  graphName = graph_Data[0]
  if (len(graph_edges) != 0 and graphNumNodes):
    print(graphNumNodes, graphName)
    all_outcome = [graphName,graphNumNodes,len(graph_edges)]
    x = testAll("ring", "SparsSuite", graph_edges , 1, graphNumNodes)
    for y in x:
      all_outcome.append(y)
    graphRes.append(all_outcome)

saveCSV(graphRes)

all_edges = read_weighted_edges_from_folder("/Users/pourdamghani/Desktop/SparseSuite")