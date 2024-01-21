#!/Library/Frameworks/Python.framework/Versions/3.10/bin/python3

import pygame
import numpy as np
import sys
import csv
import math

pygame.init()

white = (255,255,255)
black = (0,0,0)
gray = (40,40,40)
V = (128,0,255)
I = (75,0,130)
B = (0,0,205)
G = (0,128,0)
Y = (255,255,0)
O = (255,165,0)
R = (255,0,0)

(width,height) = (1400,900)
P = []

'''
P.append((20,20))
P.append((200,60))
P.append((330,255))
P.append((20,220))
P.append((20,40))
P.append((100,80))
P.append((100,120))
'''

def Pblend(P1,P2,fac):
    delta = (fac*P2[0]-fac*P1[0],fac*P2[1]-fac*P1[1])
    return (P1[0]+delta[0],P1[1]+delta[1])

def intVec(Pfloat):
    retVec = []

    if len(Pfloat) ==2:
        if type(Pfloat[0]) == float:
            return ( int(Pfloat[0]) , int(Pfloat[1]) )
    for point in Pfloat:
        #print(point)
        retVec.append((int(point[0]),int(point[1])))
    return retVec


def drawBezierShell(screen,P,t):
    Pnext = []
    lastPoint = []
    while len(P) > 2:
        Pnext = []
        for i in range(0,len(P)-1):
            Pnext.append(Pblend(P[i],P[i+1],t))
        pygame.draw.aalines(screen,white,False,intVec(Pnext),1)
        P=Pnext
    if len(Pnext) == 2:
        lastPoint = Pblend(Pnext[0],Pnext[1],t)
        #print('lastPoint: ',lastPoint)
        pygame.draw.circle(screen,G,intVec(lastPoint),5)
    return lastPoint

def calcBezier(P0,nPts):
    Pbez = []
    if len(P0) < 2:
        return []

    for j in range(0,nPts):
        fac = j/(nPts-1)
        P = P0
        while len(P) > 1:
            Pnext = []
            for i in range(0,len(P)-1):
                Pb = Pblend(P[i],P[i+1],fac)
                Pnext.append(Pb)
            P=Pnext
        Pj = Pnext[0]
        Pbez.append(Pj)
        #print(fac,len(P),Pj[0],Pj[1])

    return Pbez



print("driver: ",pygame.display.get_driver())
screen = pygame.display.set_mode((width, height),pygame.SHOWN)
pygame.display.set_caption('Gate')
screen.fill(black)
Pbez = calcBezier(P,40)


running = True
t = 0
speed = .004
delta_t = 30
last_tick = pygame.time.get_ticks()
while running:
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            running = False

    P0 = P
    t += speed
    if t > 1 or t < 0:
        speed = -1*speed

    if pygame.mouse.get_pressed()[0]:
        mousePos = pygame.mouse.get_pos()
        P.append(mousePos)

    if pygame.mouse.get_pressed()[2]:
        mousePos = pygame.mouse.get_pos()
        if len(P) >= 1:
            P.pop()


    print(t)
    screen.fill(black)
    for point in P0:
        pygame.draw.circle(screen,B,point,5)
    if len(P0) >=2:
        #print('GonnaDraw Points: ',len(P0),P0)
        pygame.draw.aalines(screen,white,False,P0,1)
        Pbez = calcBezier(P0,100*len(P0))
        pygame.draw.aalines(screen,gray,False,Pbez)
        bezPoint = drawBezierShell(screen,P0,t)
    new_tick = pygame.time.get_ticks()
    delta_tick = new_tick-last_tick
    pygame.display.flip()
    pygame.time.delay(delta_t-delta_tick)
    last_tick = new_tick




pygame.quit()
