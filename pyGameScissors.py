import myPygameVecOps as vecOp
import pygame

class pyGameScissors(pygame.sprite.Sprite):
    def __init__(self,initPos=(20,20), velocity = (0,0), radius=10, color=(0,128,0),mass=5):
        pygame.sprite.Sprite.__init__(self)
        self.image = pygame.image.load('/Users/ryanblanchard/myApplications/Python_bin/_img/scissors.png')
        self.rect = self.image.get_rect()
        self.position = initPos
