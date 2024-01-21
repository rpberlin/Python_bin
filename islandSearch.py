#!/usr/bin/env python3
# -*- coding: utf-8 -*-



grid =  [[1, 0, 1], [0, 1, 0], [1, 0, 1]]
grid = [[0, 0, 0], [0, 1, 0], [0, 0, 0]]
grid = [[1, 1, 0, 0, 1, 1],[1, 1, 0, 0, 1, 1],[1, 0, 0, 0, 1, 0],[1, 0, 0, 0, 1, 0],[1, 1, 0, 0, 1, 1],[1, 1, 0, 0, 1, 1]]
grid = [[1, 1, 0, 0, 1, 1],[1, 1, 0, 0, 1, 1],[1, 0, 0, 0, 1, 0],[1, 0, 0, 0, 1, 0],[1, 1, 0, 0, 1, 1],[1, 1, 1, 1, 1, 1]]
  
print('GRID:')
for row in grid:
    print(row)
print('\n\n')    


def cellIDfromIJ(i,j,m,n):
    return j+i*n


def getNeighborOnes(i,j,m,n,grid):
    ret = []
    if i<n-1:
        if grid[i+1][j] == 1:
            ret.append(cellIDfromIJ(i+1,j,m,n))
    if j<m-1:
        if grid[i][j+1] == 1:
            ret.append(cellIDfromIJ(i,j+1,m,n))
    if i>0:
        if grid[i-1][j] == 1:
            ret.append(cellIDfromIJ(i-1,j,m,n))
    if j>0:
        if grid[i][j-1] == 1:
            ret.append(cellIDfromIJ(i,j-1,m,n))
    return ret, len(ret)
    
    

def isInTheseListsAlready(node,graph):
    ret = []
    for i, list in enumerate(graph):
        if node in list:
            ret.append(i)
    return ret, len(ret)

def addEdge(graph, node1,node2):
    for idx, list in enumerate(graph):
        if node1 in list:
            if node2 not in list:
                graph[idx].append(node2)
    return

            
        


graph = []


n = len(grid)
m = len(grid[0])
print(n,m)
for i in range(0,n):
    for j in range(0,m):
        currNode = cellIDfromIJ(i,j,m,n)
        if grid[i][j] == 1:
            listIdxs, nLists = isInTheseListsAlready(currNode,graph)
            if nLists == 0:
                graph.append([currNode])
            else:
                for idx in listIdxs:
                    graph[idx].append(currNode)
            
        neighborOnes, nNeighborOnes = getNeighborOnes(i,j,m,n,grid)
        for neighbor in neighborOnes:
            addEdge(graph,currNode,neighbor)
        
print('GRAPH:')            
for row in graph:
    print(row)            



