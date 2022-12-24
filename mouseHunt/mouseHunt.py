#!/Library/Frameworks/Python.framework/Versions/3.10/bin/python3

import pygame
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

(width,height) = (1300,700)
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
        if pygame.mouse.get_pressed()[1]:
            mousePos = pygame.mouse.get_pos()
            newCat = mh.mhCat(mousePos)
            allCats.add(newCat)
            allSprites.add(newCat)
            catForces.append((0,0))
            print('CLICKED!: ',mousePos)
        if pygame.mouse.get_pressed()[2]:
            mousePos = pygame.mouse.get_pos()
            newMouse = mh.mhMouse(mousePos)
            allMouses.add(newMouse)
            allSprites.add(newMouse)
            mouseForces.append((0,0))
            print('CLICKED!: ',mousePos)


    mousePos = pygame.mouse.get_pos()
    theCheese.setMotion(mousePos)
    theCheese.draw(screen)
    #mousePos = (400,150)
    catForces, mouseForces = mh.calcCatAndMouseForces(allCats,allMouses,catForces,mouseForces,theCheese)
    for idx, thisCat in enumerate(allCats):
        thisCat.integrateMotion(catForces[idx],delta_t*speed,window_size)
        thisCat.draw(screen)
    for idx, thisMouse in enumerate(allMouses):
        print('MouseForce Info: ',idx,mouseForces[idx])
        thisMouse.integrateMotion(mouseForces[idx],delta_t*speed,window_size)
        thisMouse.draw(screen)
    allMouses, mouseForces = mh.catMousecollisions(allCats,allMouses,catForces,mouseForces)
    new_tick = pygame.time.get_ticks()
    delta_tick = new_tick-last_tick
    print('nSprites: ',len(allSprites),' delta_tick: ',delta_tick)
    pygame.display.flip()
    pygame.time.delay(delta_t-delta_tick)
    last_tick = new_tick
    screen.fill(black)


pygame.quit()
