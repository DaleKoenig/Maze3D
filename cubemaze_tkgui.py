#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  4 11:36:15 2019

@author: Dale Koenig
"""

import tkinter as tk
from cubemaze import rotate_mapping, direction_mapping, MazeNav3D, MazeEntity

def get_bounding_box(n,layer,i,j):
    size = 20
    spacing = 40
    base_y = (layer // 2) * (spacing + (2*n+1)*size) + i * size
    base_x = j * size
    return (base_x, base_y, base_x+size, base_y+size)

def draw_maze(maze,canv):
    canv.delete('all')
    for layer in range(1,2*n,2):
        for row in range(2*n+1):
            for col in range(2*n+1):
                val = maze[layer][row][col]
                if val == 1: # wall
                    canv.create_rectangle(*get_bounding_box(n,layer,row,col),fill='#000000')
                elif val == 2: #player
                    canv.create_oval(*get_bounding_box(n,layer,row,col),fill='#0000FF')
                elif val == 3: #goal
                    canv.create_oval(*get_bounding_box(n,layer,row,col),fill='#00FF00')
                if val != 1: #hint at paths above and below
                    bounds = get_bounding_box(n,layer,row,col)
                    if maze[layer-1][row][col] != 1: #can move up
                        canv.create_polygon([(bounds[0]+bounds[2]) // 2, bounds[1], \
                                            int(.8*bounds[0] + .2*bounds[2]), int(.7*bounds[1]+.3*bounds[3]), \
                                            int(.2*bounds[0] + .8*bounds[2]), int(.7*bounds[1]+.3*bounds[3])])
                    if maze[layer+1][row][col] != 1: #can move down
                        canv.create_polygon([(bounds[0]+bounds[2]) // 2, bounds[3], \
                                            int(.8*bounds[0] + .2*bounds[2]), int(.3*bounds[1]+.7*bounds[3]), \
                                            int(.2*bounds[0] + .8*bounds[2]), int(.3*bounds[1]+.7*bounds[3])])
def move_player(maze_game,canv,direction):
    direction = direction_mapping[direction]
    res = maze_game.attempt_move('player',direction)
    if res.name == 'goal':
        # Player won
        canv.delete('all')
        canv.create_text(50,100,text = 'YOU WIN')
    else:
        draw_maze(maze_game.maze,canv) 
                   
def rotate(maze_game,canv,rot):
    rot = rotate_mapping[rot]
    maze_game.apply_isometry(rot)
    draw_maze(maze_game.maze,canv)
                                   
root= tk.Tk()
canv = tk.Canvas(root, width=300, height=1000)
canv.grid(row=0,column=0, rowspan=6, columnspan=6)

n = 4
maze_game = MazeNav3D(n)
player = MazeEntity('player',('/\\','\\/'),(1,1,1))
goal = MazeEntity('goal',('OO','OO'), (-2,-2,-2))
maze_game.add_entity(player)
maze_game.add_entity(goal)
draw_maze(maze_game.maze,canv)

north = tk.Button(root, text='N', width=10, height=5, \
                  command = lambda : move_player(maze_game,canv,'w'))
north.grid(row=2,column=7)
west = tk.Button(root,text='W',width=10,height=5, \
                  command = lambda : move_player(maze_game,canv,'a'))
west.grid(row=3,column=6)
east = tk.Button(root,text='E',width=10,height=5, \
                  command = lambda : move_player(maze_game,canv,'d'))
east.grid(row=3,column=8)
south = tk.Button(root,text='S',width=10,height=5, \
                  command = lambda : move_player(maze_game,canv,'s'))
south.grid(row=4,column=7)
rnorth = tk.Button(root, text='RN', width=10, height=5, \
                  command = lambda : rotate(maze_game,canv,'rw'))
rnorth.grid(row=1,column=7)
rwest = tk.Button(root,text='RW',width=10,height=5, \
                  command = lambda : rotate(maze_game,canv,'ra'))
rwest.grid(row=3,column=5)
reast = tk.Button(root,text='RE',width=10,height=5, \
                  command = lambda : rotate(maze_game,canv,'rd'))
reast.grid(row=3,column=9)
rsouth = tk.Button(root,text='RS',width=10,height=5, \
                  command = lambda : rotate(maze_game,canv,'rs'))
rsouth.grid(row=5,column=7)


root.wm_attributes("-topmost", 1)
root.mainloop()
