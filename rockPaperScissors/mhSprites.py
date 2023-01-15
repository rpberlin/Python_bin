import mhVecops as vecOp
import pygame


#self.image = pygame.image.load('/Users/ryanblanchard/myApplications/Python_bin/_img/mouse-face.svg')

#self.image = pygame.image.load('/Users/ryanblanchard/myApplications/Python_bin/_img/cheese.png')

class mhSprite(pygame.sprite.Sprite):
    def __init__(self,initPos,image, velocity = (0,0), radius=10, color=(0,128,0),mass=5):
        pygame.sprite.Sprite.__init__(self)
        self.image = image
        self.delta_pos = (-1*int(self.image.get_rect()[2]/2),-1*int(self.image.get_rect()[3]/2))
        self.radius = min(self.delta_pos)
        print("### DeltaPOS: ",self.delta_pos )
        self.position = initPos
        self.velocity  = velocity
        self.color = color
        self.mass = mass
        self.tail = []

    def __str__(self):
        ret = f"Ball Position: {self.position}\n"
        ret = ret+f"Ball Velocity: {self.velocity}\n"
        ret = ret+f"Ball Radius: {self.radius}\n"
        ret = ret+f"Ball Color: {self.color}\n"
        ret = ret+f"Ball mass: {self.mass}\n"
        return ret

    def setMotion(self,position,velocity=(0,0)):
        self.position = position
        self.velocity = velocity

    def integrateMotion(self, Force, delta_t, window_size, rCoeffNorm = .6, rCoeffTang = .6):

        if self.position[0] < 0:
            self.position = (0,self.position[1])
            self.velocity = (-1*rCoeffNorm*self.velocity[0],rCoeffTang*self.velocity[1])
        if self.position[1] < 0:
            self.position = (self.position[0],0)
            self.velocity = (rCoeffTang*self.velocity[0],-1*rCoeffNorm*self.velocity[1])
        if self.position[0] > window_size[0]:
            self.position = (window_size[0],self.position[1])
            self.velocity = (-1*rCoeffNorm*self.velocity[0],rCoeffTang*self.velocity[1])
        if self.position[1] > window_size[1]:
            self.position = (self.position[0],window_size[1])
            self.velocity = (rCoeffTang*self.velocity[0],-1*rCoeffNorm*self.velocity[1])


        self.velocity = vecOp.floatVectorAdd( self.velocity , vecOp.floatVectorScale(Force,delta_t/self.mass) )
        self.position = vecOp.floatVectorAdd( self.position , vecOp.floatVectorScale(self.velocity,delta_t) )
        self.tail.append(vecOp.intVector(self.position))
        if len(self.tail) > 60:
            del self.tail[0]
        return

    def removeTail(self):
        self.tail = []
        return


    def draw(self,screen):
        #pygame.draw.circle(screen,self.color,vecOp.intVector(self.position),self.radius)
        if len(self.tail) > 2:
            pygame.draw.aalines(screen,self.color,False,self.tail,3)
        screen.blit(self.image,vecOp.intVectorAdd(self.position,self.delta_pos))
        pass

class mhMouse(mhSprite):
    def __init__(self,initPos):
        image = pygame.image.load('/Users/ryanblanchard/myApplications/Python_bin/mouseHunt/_img/mouse-face.svg')
        mhSprite.__init__(self,initPos,image)



class mhCat(mhSprite):
    def __init__(self,initPos):
        image = pygame.image.load('/Users/ryanblanchard/myApplications/Python_bin/mouseHunt/_img/cat-face.png')
        mhSprite.__init__(self,initPos,image)

class mhCheese(mhSprite):
    def __init__(self,initPos):
        image = pygame.image.load('/Users/ryanblanchard/myApplications/Python_bin/mouseHunt/_img/cheese.png')
        mhSprite.__init__(self,initPos,image)

catSameG = .04  #how hard cats run away from other cats
mouseSameG = .04 #how hard mice run away from each mice
catMouseG = .05  #how hard the cat run towards the mouse
mouseCatG  = .08 #how hard the mouse runs away from the cats
cheeseGravity = .00003 #how strong the spring is between the mouse and the cheese
def calcCatAndMouseForces(cats,mouses,catForces,mouseForces,theCheese):
    for i, cat_i in enumerate(cats):
        force_cat = (0,0)
        for j, cat_j in enumerate(cats):
            if i != j:
                force_ij, d2_ij = vecOp.invDistWithComponents(cat_j.position,cat_i.position,catSameG)
                force_cat = vecOp.floatVectorAdd(force_cat,force_ij)
        for k, mouse_k in enumerate(mouses):
                force_ik, d2_ik = vecOp.invDistWithComponents(cat_i.position,mouse_k.position,catMouseG)
                force_cat = vecOp.floatVectorAdd(force_cat,force_ik)
        catForces[i] = force_cat

    for i, mouse_i in enumerate(mouses):
        force_mouse = (0,0)
        force_cheese = vecOp.floatVectorDiff(mouse_i.position,theCheese.position,cheeseGravity)
        force_mouse = vecOp.floatVectorAdd(force_mouse,force_cheese)
        for k, cat_k in enumerate(cats):
                force_ik, d2_ik = vecOp.invDistWithComponents(cat_k.position,mouse_i.position,mouseCatG)
                force_mouse = vecOp.floatVectorAdd(force_mouse,force_ik)
        for j, mouse_j in enumerate(mouses):
            if i != j:
                force_ij, d2_ij = vecOp.invDistWithComponents(mouse_j.position,mouse_i.position,mouseSameG)
                force_mouse = vecOp.floatVectorAdd(force_mouse,force_ij)
        mouseForces[i] = force_mouse
    return catForces, mouseForces

def catMousecollisions(cats,mouses,catForces,mouseForces):
    mouseKillList = []
    for i, mouse_i in enumerate(mouses):
        for k, cat_k in enumerate(cats):
                 dist2 = vecOp.intDist2(mouse_i.position,cat_k.position)
                 if dist2 < cat_k.radius*mouse_i.radius:
                     mouseKillList.append(mouse_i)
                     break
    for mouse in mouseKillList:
        mouses.remove(mouse)
        mouseForces.pop()

    return mouses,mouseForces
