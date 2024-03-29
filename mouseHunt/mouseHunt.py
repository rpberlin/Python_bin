#!/Library/Frameworks/Python.framework/Versions/3.10/bin/python3

import pygame
from pygame.locals import *
import mhSprites as mh
import mhVecops as vecOp
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

(width,height) = (1600,900)
window_size = (width,height)
P = []

P.append((20,20))
P.append((200,60))
P.append((330,255))
P.append((20,220))
P.append((20,40))
P.append((100,80))
P.append((100,120))


screen = pygame.display.set_mode((width, height),flags=pygame.SHOWN)
pygame.display.set_caption('MouseHunter - 1 Ball')
screen.fill(black)
running = True



allSprites = pygame.sprite.Group()
allMouses  = pygame.sprite.Group()
allCats    = pygame.sprite.Group()
allDogs    = pygame.sprite.Group()
firstMouse = mh.mhMouse((500,130))
theCheese = mh.mhCheese((500,130))
allMouses.add(firstMouse)
allSprites.add(firstMouse)
fps = 30
speed = .3
delta_t = int(1000/fps)
last_tick = pygame.time.get_ticks()
mouseForces = []
catForces = []
for _ in allMouses:
    mouseForces.append((0,0))
for _ in allCats:
    catForces.append((0,0))


while running:
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            running = False
        if pygame.mouse.get_pressed()[0]:
            mousePos = pygame.mouse.get_pos()
            firstMouse.setMotion(mousePos,(0,0))
            firstMouse.removeTail()
        #if pygame.mouse.get_pressed()[1]:
        if pygame.key.get_pressed()[K_c]:
            mousePos = pygame.mouse.get_pos()
            newCat = mh.mhCat(mousePos)
            allCats.add(newCat)
            allSprites.add(newCat)
            catForces.append((0,0))
            print('CLICKED!: ',mousePos)
        #if pygame.mouse.get_pressed()[2]:
        if pygame.key.get_pressed()[K_m]:
            mousePos = pygame.mouse.get_pos()
            newMouse = mh.mhMouse(mousePos)
            allMouses.add(newMouse)
            allSprites.add(newMouse)
            mouseForces.append((0,0))
            print('CLICKED!: ',mousePos)

        if pygame.key.get_pressed()[K_d]:
            mousePos = pygame.mouse.get_pos()
            newDog = mh.mhDog(mousePos)
            allDogs.add(newDog)
            allSprites.add(newDog)
            print('CLICKED!: ',mousePos)



    mousePos = pygame.mouse.get_pos()
    theCheese.setMotion(mousePos)
    theCheese.draw(screen)
    #mousePos = (400,150)
    mh.calcCatAndMouseForces(allCats,allMouses,allDogs,theCheese)
    for idx, thisCat in enumerate(allCats):
        thisCat.integrateMotion(delta_t*speed,window_size)
        thisCat.draw(screen)
    for idx, thisMouse in enumerate(allMouses):
        thisMouse.integrateMotion(delta_t*speed,window_size)
        thisMouse.draw(screen)
    for idx, thisDog in enumerate(allDogs):
        thisDog.integrateMotion(delta_t*speed,window_size)
        thisDog.draw(screen)
    allMouses, allCats, allDogs = mh.catMousecollisions(allMouses,allCats,allDogs, theCheese)
    new_tick = pygame.time.get_ticks()
    delta_tick = new_tick-last_tick
    print('nSprites: ',len(allSprites),'nMouses: ',len(allMouses),'nCats: ',len(allCats),'nDogs: ',len(allDogs),' delta_tick: ',delta_tick)
    pygame.display.flip()
    pygame.time.delay(max(0,delta_t-delta_tick))
    last_tick = new_tick
    screen.fill(black)


pygame.quit()
