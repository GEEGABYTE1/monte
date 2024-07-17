import pygame 
import sys 

from time import process_time, time 
from rocketpy import Environment, SolidMotor, Rocket, Flight
import numpy as np 
from numpy.random import normal, choice
from IPython.display import display
import matplotlib as mpl 
import matplotlib.pyplot as plt

pygame.init()  

screen_width = 800
screen_height = 600 
screen = pygame.display.set_mode((screen_width, screen_height))
pygame.display.set_caption('RocketPY GUI')

running = True 

#colors
WHITE = (255, 255, 255)
BLACK = (0, 0, 0)

#fonts
title_font = pygame.font.Font(None, 50)
subtitle_font = pygame.font.Font(None, 32)
text_font = pygame.font.Font(None, 18)

#rendering init
title_surface = title_font.render('Rocketpy GUI', True, BLACK)
subtitle_surface = subtitle_font.render('monte carlo sims done easy', True, BLACK)
text_surface = text_font.render('Some additional text.', True, BLACK)


#rect and pos
title_rect = title_surface.get_rect(center=(screen_width // 2, 100))
subtitle_rect = subtitle_surface.get_rect(center=(screen_width // 2, 150))
text_rect = text_surface.get_rect(center=(screen_width // 2, 300))

while running:

    #setting matplotlib stuff 
    mpl.rcParams["figure.figsize"] = [8, 5]
    mpl.rcParams["figure.dpi"] = 120 
    mpl.rcParams["font.size"] = 14 
    mpl.rcParams["legend.fontsize"] = 14 
    mpl.rcParams["figure.titlesize"] = 14




    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            running = False 
        
    
    screen.fill(WHITE)

    #content 
    screen.blit(title_surface, title_rect)
    screen.blit(subtitle_surface, subtitle_rect)
    screen.blit(text_surface, text_rect)


    #Updating 
    pygame.display.flip()

pygame.quit()
sys.exit()