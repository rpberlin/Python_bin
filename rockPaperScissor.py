#!/Library/Frameworks/Python.framework/Versions/3.10/bin/python3

import pygame
import myPygameBall as pb
import myPygameVecOps as vecOp
import pyGameScissors as scis
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


(width,height) = (900,500)
window_size = (width,height)
screen = pygame.display.set_mode((width, height),flags=pygame.SHOWN)


newScis = scis.pyGameScissors()
delta_t =1000

running = True
posx = 40
posy = 40
t = 0
speed = .002
while running:
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            running = False

    screen.fill(gray)
    t+=1
    position = (posx+t, posy+t)
    newScis.position = position
    screen.blit(newScis.image,newScis.position)
    pygame.display.flip()
    pygame.time.delay(delta_t)



pygame.quit()
