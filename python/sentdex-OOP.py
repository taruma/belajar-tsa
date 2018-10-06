# -*- coding: utf-8 -*-
"""
Created on Fri Sep 28 08:41:49 2018

@author: tarum
"""
import random
import pygame
from blob import Blob

STARTING_BLUE_BLOBS = 10
STARTING_RED_BLOBS = 3

WIDTH = 800
HEIGHT = 600
WHITE = (255, 255, 255)
BLUE = (0, 0, 255)
RED = (255, 0, 0)

game_display = pygame.display.set_mode((WIDTH, HEIGHT))
pygame.display.set_caption("BLOB WORLD")
clock = pygame.time.Clock()

class BlueBlob(Blob): #BLOB = parent, super, base-class
    def __init__(self, color, x_boundary, y_boundary):
        super().__init__(color, x_boundary, y_boundary)
#        Blob.__init__(self, color, x_boundary, y_boundary)
        self.color = BLUE
    
    def move_fast(self):
        self.x += random.randrange(-7,7)
        self.y += random.randrange(-7,7)


def draw_environment(blob_list):
    game_display.fill(WHITE)
    
    for blob_dict in blob_list:
        for blob_id in blob_dict:
            blob = blob_dict[blob_id]
            pygame.draw.circle(game_display, blob.color, [blob.x, blob.y], blob.size)
            blob.move_fast()
            blob.check_bounds()
#            if blob.x < 0: blob.x = 0
#            elif blob.x > blob.x_boundary: blob.x = blob.x_boundary
#            
#            if blob.y < 0: blob.y = 0
#            elif blob.y > blob.y_boundary: blob.y = blob.y_boundary
    pygame.display.update()
    
    
def main():
#    red_blob = Blob(RED)
    blue_blobs = dict(enumerate([BlueBlob(BLUE, WIDTH, HEIGHT) for i in range(STARTING_BLUE_BLOBS)]))
    red_blobs = dict(enumerate([BlueBlob(RED, WIDTH, HEIGHT) for i in range(STARTING_RED_BLOBS)]))
#    print(blue_blobs)
    while True:
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                pygame.quit()
#                quit()
        
        draw_environment([blue_blobs, red_blobs])
        clock.tick(60)
#        print(red_blob.x, red_blob.y)

if __name__ == '__main__':
    main()