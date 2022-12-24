#!/Library/Frameworks/Python.framework/Versions/3.10/bin/python3

import pygame
import numpy as np
import sys
import csv
import math

pygame.init()

white = (255,255,255)
black = (0,0,0)
V = (128,0,255)
I = (75,0,130)
B = (0,0,205)
G = (0,128,0)
Y = (255,255,0)
O = (255,165,0)
R = (255,0,0)

(width,height) = (400,265)
P = []

P.append((20,20))
P.append((200,60))
P.append((330,255))
P.append((20,220))
P.append((20,40))
P.append((100,80))
P.append((100,120))


def Pblend(P1,P2,fac):
    delta = (int(fac*(P2[0]-P1[0])),int(fac*(P2[1]-P1[1])))
    return (P1[0]+delta[0],P1[1]+delta[1])


def drawMultiLine(screen,P):
    pygame.draw.aalines(screen,white,False,P,1)
    #for i in range(0,len(P)-1):
    #    pygame.draw.aaline(screen,white,P[i],P[i+1],1)
    return

def drawBezierShell(screen,P,t):
    while len(P) > 2:
        Pnext = []
        for i in range(0,len(P)-1):
            Pnext.append(Pblend(P[i],P[i+1],t))
        drawMultiLine(screen,Pnext)
        P=Pnext
    return

def calcBezier(P0,nPts):
    Pbez = []
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
        print(fac,len(P),Pj[0],Pj[1])

    return Pbez



print("driver: ",pygame.display.get_driver())
screen = pygame.display.set_mode((width, height),pygame.SHOWN)
pygame.display.set_caption('Gate')
screen.fill(black)
Pbez = calcBezier(P,40)


running = True
t = 0
speed = .002

while running:
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            running = False
    pygame.display.flip()
    P0 = P
    t += speed
    if t > 1 or t < 0:
        speed = -1*speed

    print(t)
    screen.fill(black)
    drawMultiLine(screen,P0)
    drawMultiLine(screen,Pbez)
    drawBezierShell(screen,P0,t)
    pygame.time.delay(30)


pygame.quit()
