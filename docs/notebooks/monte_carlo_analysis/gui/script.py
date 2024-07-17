import pygame
import sys
import matplotlib as mpl

pygame.init()

screen_width = 1000
screen_height = 1000
screen = pygame.display.set_mode((screen_width, screen_height))
pygame.display.set_caption('RocketPY GUI')

running = True

# Colors
WHITE = (255, 255, 255)
BLACK = (0, 0, 0)

# Image
image_path = './images/rocket.jpg'
image = pygame.image.load(image_path)

# Desired image size
desired_width = 400
aspect_ratio = image.get_height() / image.get_width()
scaled_image = pygame.transform.scale(image, (desired_width, int(desired_width * aspect_ratio)))
image_rect = scaled_image.get_rect(center=(screen_width // 2, screen_height // 2))

button_text = 'Get Started'

# Fonts
title_font = pygame.font.Font(None, 50)
subtitle_font = pygame.font.Font(None, 32)
text_font = pygame.font.Font(None, 18)
button_font = pygame.font.Font(None, 36)

# Rendering
title_surface = title_font.render('Rocketpy GUI', True, BLACK)
subtitle_surface = subtitle_font.render('Monte Carlo sims made easy', True, BLACK)
text_surface = text_font.render('Made by Jaival Patel @ 2024', True, BLACK)
button_text_surface = button_font.render(button_text, True, WHITE)



title_rect = title_surface.get_rect(center=(screen_width // 2, 100))
subtitle_rect = subtitle_surface.get_rect(center=(screen_width // 2, 150))
text_rect = text_surface.get_rect(center=(screen_width // 2, 900))
button_rect = pygame.Rect(400, 400, 200, 50)
button_text_rect = button_text_surface.get_rect(center=button_rect.center)

button_color = (0, 128, 0)


button_clicked = False



while running:
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            running = False
        elif event.type == pygame.MOUSEBUTTONDOWN:
            if event.button == 1: #checking left mouse button clicked
                if button_rect.collidepoint(event.pos):
                    button_clicked = not button_clicked

    screen.fill(WHITE)


    if button_clicked:
        button_text_surface = button_font.render('Clicked!', True, WHITE)
    else:      
        screen.blit(scaled_image, image_rect)
        screen.blit(title_surface, title_rect)
        screen.blit(subtitle_surface, subtitle_rect)
        screen.blit(text_surface, text_rect)
        pygame.draw.rect(screen, button_color, button_rect)
        screen.blit(button_text_surface, button_text_rect)
        button_text_surface = button_font.render('Get Started!', True, WHITE)

    button_text_rect = button_text_surface.get_rect(center=button_rect.center)

    pygame.display.flip()

pygame.quit()
sys.exit()
