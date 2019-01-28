#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Cubic maze game.

The goal of this game is to create a 3d maze to practice and showcase python
programming abilities.

I also put a bit more detail than needed into the class used for rotating cubes
(and thus used to rotate the maze)

Created on Fri Jan 25 16:47:20 2019

@author: Dale Koenig
"""
import re
import numpy as np

class CubeIsometry:
    """A class for representing isometries of the cube, with operators
    overloaded to allow easy multiplication, powers, and application to nxnxn
    array like objects.
    
    During initialization, the isometry may be specified by a string giving a
    permutation of the three axes 'i', 'j', and 'k'.  So for example 'i -k j' 
    is parsed as
    i -> i
    j -> -k
    k -> j
    and therefore corresponds to rotating the k axis 90 degrees to the j axis.
    
    Properties:
        isom -- A string representation of the isometry.  Used to initialize.
        signs --  A list specifying whether each of i, j, k gets sent to the
            positive or negative target vector in the permutation.  1 if
            positive, -1 if negative.
        permutation -- A list specifying the permutation applied to the i, j, 
        and k axes.
    """
    

    def __init__(self, isom):
        """
        Args:
            isom -- A string representing the isometry to initialize to.  See
                the class docstring for guidelines on how these strings are
                defined.
        """
        self.isom = isom
        # Lets sanitize isom with a regex check followed by checking if all of
        # i, j, k  are contained in it.
        sanitycheck = re.compile('-?[ijk] -?[ijk] -?[ijk]$')
        if not sanitycheck.match(isom):
            raise ValueError('Invalid isometry specification "{}"'.format(isom)) 
        if not ('i' in isom and 'j' in isom and 'k' in isom):
            raise ValueError('Invalid isometry specification "{}"'.format(isom)) 
        axes = isom.split()
        signs = [1, 1, 1]
        permutation = [-1,-1,-1]
        for index, axis in enumerate(axes):
            if axis[0] == '-':
                signs[index] = -1
            if axis[-1] == 'i':
                permutation[index] = 0
            elif axis[-1] == 'j':
                permutation[index] = 1
            elif axis[-1] == 'k':
                permutation[index] = 2
        self.signs = signs
        self.permutation = permutation
        
    def __mul__(self, other):
        # Calculate signs of the product by multiplying signs of the two objects
        prod_signs = [self.signs[other.permutation[i]]*other.signs[i] for i in range(3)]
        # Calculate the permutation of the product by chaining
        # See 'one line representation' of permutations and how to multiply them
        prod_permutation = [self.permutation[other.permutation[i]] for i in range(3)]
        # Now begin to compose the final string.
        prefixes = ['-' if sign == -1 else '' for sign in prod_signs]
        postfixes = ['ijk'[target] for target in prod_permutation]
        result = ' '.join([prefixes[i] + postfixes[i] for i in range(3)])
        return CubeIsometry(result)
    
    def __eq__(self, other):
        return(self.isom == other.isom)
        
    def __str__(self):
        return self.isom
    
    def __pow__(self,exponent):
        """We calulate the power in the naive way.  This cannot get too 
        inefficient because the largest possible exponent is 12.  Indeed, the
        possible orders of elements in the group are 1,2,3,4, and 6, and the lcm
        of these orders is 12.
        
        This can be used to calculate the inverse by using -1 as the exponent.
        """
        exponent = exponent % 12 # 12 is the lcm of all possible orders in the group
        base = CubeIsometry('i j k')
        for i in range(exponent):
            base = self * base
        return base
    
    def get_order(self):
        """Calculates and returns the order of the isometry, that is the
         minimum positive n such that x^n == CubeIsometry('i j k')
         """
        for i in range(1,25):
            if (self ** i) == CubeIsometry('i j k') :
                return i
    
    def preserves_orientation(self):
        """Calculates and returns True if the isometry preserves the
        orientation of the cube, and False otherwise.
        
        Each -1 in self.signs represents a reflection, and self.permutation 
        contributes a -1 iff the permutation is odd.
        """
        sign_prod = self.signs[0] * self.signs[1] * self.signs[2]
        perm_sign = (-1 if self.permutation in [[0,2,1],[1,0,2],[2,1,0]] else 1)
        if sign_prod * perm_sign == 1:
            return True
        else:
            return False
        
    def act_on_cube(self,cube):
        """Apply the isometry to a cube (an n1 x n2 x n3 numpy array) and return
        the resulting array.
        """
        # First figure out which axes are reflected by the isometry and flip those
        flip_axes = []
        if self.signs[0] == -1:
            flip_axes.append(0)
        if self.signs[1] == -1:
            flip_axes.append(1)
        if self.signs[2] == -1:
            flip_axes.append(2)
        cube = np.flip(cube,flip_axes)
        
        # Now permute according to the isometry permutation
        cube = np.transpose(cube,axes=self.permutation)
        
        return cube
    
    def act_on_point(self,point):
        """Apply the isometry to a point lying in a cube."""
        i,j,k = point
        if self.signs[0] == -1:
            i = -1 -i
        if self.signs[1] == -1:
            j = -1 -j
        if self.signs[2] == -1:
            k = -1 -k
        perm_inverse = (self.permutation.index(0), \
                        self.permutation.index(1), \
                        self.permutation.index(2))
        unpermuted_coords = (i,j,k)
        permuted_coords = (unpermuted_coords[perm_inverse[0]], \
                           unpermuted_coords[perm_inverse[1]], \
                           unpermuted_coords[perm_inverse[2]])
        return permuted_coords
        
 
def get_all_isometries():
    """Return a list containing every possible CubeIsometry object.  
    Useful for testing.
    """
    isoms = []
    possible_signs = [[i,j,k] for i in [-1,1] for j in [-1,1] for k in [-1,1]]
    possible_perms = [[0,1,2],[0,2,1],[1,0,2],[1,2,0],[2,0,1],[2,1,0]]
    for signs in possible_signs:
        for perm in possible_perms:
            prefixes = ['-' if sign == -1 else '' for sign in signs]
            postfixes = ['ijk'[target] for target in perm]
            result = ' '.join([prefixes[i] + postfixes[i] for i in range(3)])
            isoms.append(CubeIsometry(result))
    return isoms
    
def cube_isometry_test():
    """Run code validation checks on CubeIsometry.
    """
    # Default initialization test
    x = CubeIsometry()
    if(x.signs != [1,1,1]):
        print('Initialized identity but got signs {}'.format(x.signs))
    if(x.permutation != [0,1,2]):
        print('Initialized identity but got permutation {}'.format(x.permutation))
    
    # Initialize with some negatives
    x = CubeIsometry('i -k -j')
    if x.signs != [1, -1, -1]:
        print('Initialized CubeIsometry("i -k -j") but got signs {}'.format(x.signs))
    if x.permutation != [0,2,1]:
        print('Initialized CubeIsometry("i -k -j") but got permutation {}'.format(x.permutation))
        
    # Attempt invalid isometry
    failed = False
    try:
        x = CubeIsometry('j j k')
    except (ValueError):
        failed = True
    finally:
        if not failed:
            print('Expected initialiation of isometry "j j k" to fail but it succeeded')
            
    # Attempt invalid input
    failed = False
    try:
        x = CubeIsometry('i -j -ikji')
    except (ValueError):
        failed = True
    finally:
        if not failed:
            print('Expected initialiation of isometry "i -j -iiji" to fail but it succeeded')
    
    # Attempt some products.
    if not (CubeIsometry('k j i') == CubeIsometry('i k j') * CubeIsometry('j k i')):
        print('Product CubeIsometry("i k j") * CubeIsometry("j k i") gave unexpected result:')
        print(CubeIsometry('i k j') * CubeIsometry('j k i'))
        
    if not (CubeIsometry('j k -i') == CubeIsometry('-i -k -j') * CubeIsometry('-k -j i')):
        print('Product CubeIsometry("-i -k -j") * CubeIsometry("-k -j i") gave unexpected result:')
        print(CubeIsometry('-i -k -j') * CubeIsometry('-k -j i'))
    
    # Test order calculation and get_all_isometries()  
    correctorders = [1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3, \
                     3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,6,6,6,6,6,6,6,6]
    if (correctorders != sorted(x.get_order() for x in get_all_isometries())):
        print('Error in order calculation or in get_all_isometries()')
        
    # Test cube action
    arr = np.array([[[1,2,3],[4,5,6],[7,8,9]],          \
                    [[10,11,12],[13,14,15],[16,17,18]], \
                    [[19,20,21],[22,23,24],[25,26,27]]])
    if not (arr[::-1,:,:] == CubeIsometry('-i j k').act_on_cube(arr)).all():
        print('Error in cube isometry acting by reflection')
    transpose_test = True
    for i in range(3):
        for j in range(3):
            for k in range(3):
                if arr[j,k,i] != CubeIsometry('k i j').act_on_cube(arr)[i,j,k]:
                    transpose_test = False
    if not transpose_test:
        print('Fail in cube isometry acting by transposition')
        
    # Test point action
    if CubeIsometry('k -i j').act_on_point((1,2,3)) != (-2,3,1):
        print('Error in action of CubeIsometry("k -i j") on point (1,2,3)')
    print('Testing complete.')

DISPLAY_OPTIONS_MED = {'layers_per_row' : 7, 
                       'wall_symbol' : ('##','##'),
                        'empty_symbol' : ('  ','  '),
                        'separator' : '   '}
DISPLAY_OPTIONS_DENSE = {'layers_per_row' : 7, 
                         'wall_symbol' : ('#',),
                         'empty_symbol' : (' ',),
                         'separator' : '  '}

DISPLAY_OPTIONS_LARGE = {'layers_per_row' : 5, 
                         'wall_symbol' : ('###','###','###'),
                         'empty_symbol' : ('   ','   ','   '),
                         'separator' : '     '}
        
class Maze3D():
    """Create a 3d maze stored as a numpy array.  The maze can be rotated by 
    a CubeIsometry object.
    
    An n x n x n maze is stored as a 2n+1 x 2n+1 + 2n+1 numpy array.  Coordinates
    with all of i, j, k odd represent the maze itself, while coordinates with
    one of i, j, k even represent the walls.  
    
    Properties:
        maze -- a numpy uint8 array storing the maze itself, with 0 for walkable
                    area and 1 for walls                    
        n -- the length of each walkable dimension.  The maze itself, including
                 walls, is 2n+1 in each dimension
        seed -- the seed used to generate the maze.  Generated randomly if given
                    as None
    """
    
    def __init__(self, n, openness = 0.0, seed = None, display_options='dense'):
        """
        Args:
            n -- the size of each dimension not including walls
            
            openness -- a float between 0.0 and 1.0 indicating how restricted the
                maze is.  A value of 0.0 means there is exactly one path connecting
                each pair of points, while a value of 1.0 means there are no 
                walls obstructing movement, and only exterior walls and "pillars"
                at each i,j,k where all of i,j,k are even.  
                
                The maze is first created with unique
                paths connecting pairs of points, and then each wall is removed 
                with probability 'openness'
                
            seed -- a random seed to be passed to numpy
            
            display_options -- a dictionary can be passed to specify custom options.
                preset options can be passed with 'med', 'dense', or 'large'           
        """
        self.maze = np.ones(((2*n+1),(2*n+1),(2*n+1)),dtype=np.uint8)
        self.n = n
        self.seed = seed
        
        if display_options == 'med':
            self.display_options = DISPLAY_OPTIONS_MED
        elif display_options == 'dense':
            self.display_options = DISPLAY_OPTIONS_DENSE
        elif display_options == 'large':
            self.display_options = DISPLAY_OPTIONS_LARGE
        else:
            self.display_options = display_options
            
        np.random.seed(seed=seed)
        
        # Iterate through all odd coordinates and store them in dictionary
        # We use this dictionary to keep track of which maze coordinates have
        # been connected by paths to the start position (which coordinates are
        # 'in the maze').
        base_coords = {}
        for i in range(1,2*n,2):
            for j in range(1,2*n,2):
                for k in range(1,2*n,2):
                    self.maze[i,j,k] = 0
                    # (i,j,k) has not yet been connected to the start by a path
                    # so we set its dictionary entry to 0
                    base_coords[(i,j,k)] = 0
        # Choose (1,1,1) as the start position and set its dicitonary entry to 1
        base_coords[(1,1,1)] = 1
        # Keep track of how many spots are left to connect to the start position
        not_added_count = n ** 3 - 1
        possible_additions = set(self.get_neighbors((1,1,1)))
        for _ in range(not_added_count):
            # Choose the next coordinate to connect to the maze
            rand_pos = np.random.choice(len(possible_additions))
            next_addition = list(possible_additions)[rand_pos]
            possible_additions.remove(next_addition)
            
            # Find the possible connections, requiring them to already be in the maze
            addition_nbrs = [x for x in self.get_neighbors(next_addition) if base_coords[x] == 1]
            connection_coord = addition_nbrs[np.random.choice(len(addition_nbrs))]
            
            # At this point 'connection_coord' is a coordinate in the maze, with a
            # neighbor 'next_addition' that is not yet connected to the maze
            # We now mark next_addition as being in the maze, and open up the wall
            # connecting it to connection_coord
            base_coords[next_addition] = 1
            wall_i, wall_j, wall_k = self.get_wall(next_addition,connection_coord)
            self.maze[wall_i, wall_j, wall_k] = 0
            
            # Now we update 'possible_additions'
            for new_nbr in self.get_neighbors(next_addition):
                if base_coords[new_nbr] == 0:
                    possible_additions.add(new_nbr)
        
        # Next we iterate through all interior walls, removing each with a 
        # probability 'openness'
        # We iterate over everything in the interval [1,2n] x [1,2n] x [1,2n]
        # to simplify code, which takes longer but does not add complexity.
        # Skip coordinates in which two or more coordinates are even, 
        # since these do not obstruct movement.
        if openness > 0.0:
            for i in range(1,2*n):
                for j in range(1,2*n):
                    for k in range(1,2*n):
                        # confirm at least two coordinates are odd.
                        if ((i%2) + (j%2) + (k%2) >= 2) and np.random.uniform() < openness:
                            self.maze[i,j,k] = 0
                            
    def convert_pos(self, pos):
        """Convert an (x,y,z) position tuple to a tuple of strings"""
        x,y,z = pos
        if self.maze[x,y,z] == 0:
            return self.display_options['empty_symbol']
        elif self.maze[x,y,z] == 1:
            return self.display_options['wall_symbol']

    def convert_row(self,i,j):
        """Convert a row of the maze to a list of strings.  Joining these strings
        by \n characters will display the row."""
        display_pieces = [self.convert_pos((i,j,k)) for k in range(len(self.maze[i,j]))]
        display_rows = len(display_pieces[0])
        # The second line below iterates through the rows of text to display
        # corresponding to the single maze row (each maze row is displayed as
        # multiple lines of ascii characters)
        # The first joins together all the text pieces making up that text row
        text_rows = [''.join(piece[text_row_index] for piece in display_pieces) \
                     for text_row_index in range(display_rows)]
        return text_rows
        
    
    def display_layer_range(self,start_i,end_i):
        """Display the layers only from 'start_i' to 'end_i' all on one line"""
        m = len(self.display_options['self_symbol'])
        for j in range(len(self.maze[0])):
                text_entries = [self.convert_row(i,j) for i in range(start_i,end_i)]
                # We must switch the inner and outer index of text so that
                # the layer is determined by the inner index, and the row of text
                # by the outer index
                text_entries = [[text_entries[i][j] for i in range(len(text_entries))] 
                                    for j in range(m)]
                display_rows = [self.display_options['separator'].join(row) for row in text_entries]
                curr_display ='\n'.join(display_rows)
                print(curr_display)
        
    def display_all(self):
        """Display all layers in a grid, with self.display_options['layers_per_row'] 
        determining how many layers to display side by side.  Note that if 
        self.display_options['layers_per_row'] is too large the output may be 
        incomprehensible.
        """
        m = len(self.display_options['wall_symbol'])
        for start_i in range(0,len(self.maze), self.display_options['layers_per_row']):
            end_i = min(start_i + self.display_options['layers_per_row'],len(self.maze))
            for j in range(len(self.maze[0])):
                text_entries = [self.convert_row(i,j) for i in range(start_i,end_i)]
                # We must switch the inner and outer index of text so that
                # the layer is determined by the inner index, and the row of text
                # by the outer index
                text_entries = [[text_entries[i][j] for i in range(len(text_entries))] 
                                    for j in range(m)]
                display_rows = ['     '.join(row) for row in text_entries]
                curr_display ='\n'.join(display_rows)
                print(curr_display)
            print('\n\n')
                           
    def get_neighbors(self,coords):
        """Given a coordinate tuple 'coords', returns a list of neighbors
        of coords in the maze.
        """
        i,j,k = coords
        nbrs = []
        if i >= 3:
            nbrs.append((i-2,j,k))
        if i <= 2 * self.n - 3:
            nbrs.append((i+2,j,k))
        if j >= 3:
            nbrs.append((i,j-2,k))
        if j <= 2 * self.n - 3:
            nbrs.append((i,j+2,k))
        if k >= 3:
            nbrs.append((i,j,k-2))
        if k <= 2 * self.n - 3:
            nbrs.append((i,j,k+2))
        return nbrs
            
    def get_wall(self, coord1, coord2):
        """Given coordinates of two neighboring maze positions, finds the 
        coordinates of the wall between them.
        """
        i1,j1,k1 = coord1
        i2,j2,k2 = coord2
        res = ((i1+i2)//2, (j1+j2)//2, (k1+k2)//2)
        return res


#'self_symbol' : ('---','|!|','---'),
        
class MazeEntity():
    """An entity within the maze that can possibly move around or be interacted \
    with.
    """
    
    def __init__(self, name, display, coords):
        """
        Args:
            name -- the name of the object.
            display -- a tuple of strings describing how to display the object
                this should be the same length and size as the values for
                display_options['wall_symbol'] and display_options['empty_symbol']
                in the maze game object.
            coords -- the initial location of the object in the maze.
        """
        self.name = name
        self.display = display
        self.coords = coords

class MazeNav3D(Maze3D):
    """A Maze with entities that can be moved around but cannot pass through 
    walls or through each other.
    """
    def __init__(self, 
                 n, 
                 openness = 0.0, 
                 seed = None,
                 display_options = 'dense'
                 ):
        """
        Args:
            n -- the size of each dimension not including walls   
            
            openness -- see Maze3D
            
            seed -- a random seed to be passed to numpy
            
            display_options -- a dictionary can be passed to specify custom options.
                preset options can be passed with 'med', 'dense', or 'large' 
        """
        Maze3D.__init__(self,n,openness,seed, display_options)
        # Making this a dictionary for easy modification later
        # This will control how to display the maze in ascii
        self.next_entity = 2
        self.entities = {}
    
    def apply_isometry(self, isom):
        """Applies a CubeIsometry to the maze and to each entity.
        """
        self.maze = isom.act_on_cube(self.maze)
        for item in self.entities:
            self.entities[item].coords = isom.act_on_point(self.entities[item].coords)
    def get_entity(self,name):
        """Find an entity by name"""
        for key,entity in self.entities.items():
            if entity.name == name:
                return entity
            
    def add_entity(self, ent):
        """Adds a new entity to the maze (for example, a player, or a 
        monster) that can be moved by name or interacted with.  Duplicate entity
        names are forbidden.  Overwrites any wall entity that exists at the same
        location.
        
        Args:
            ent -- A MazeEntity object specifying the entity to add.
        """
        if ent.name in self.entities:
            raise ValueError('Object {} could not be created in the maze since \
                             an object with this name already exists'.format(ent.name))
        self.entities[self.next_entity] = ent
        i,j,k = ent.coords
        self.maze[i,j,k] = self.next_entity
        self.next_entity += 1
        
    def convert_pos(self, pos):
        """Convert an (x,y,z) position tuple to a tuple of strings"""
        x,y,z = pos
        if self.maze[x,y,z] == 0:
            return self.display_options['empty_symbol']
        elif self.maze[x,y,z] == 1:
            return self.display_options['wall_symbol']
        else:
            ent = self.entities[self.maze[x,y,z]]
            return ent.display
    
    
    def attempt_move(self, ent_name, direction):
        """Moves the entity in the direction given, if there is no wall blocking
        the way.  Returns the name of the entity preventing movement if any (
        'wall' if it is a wall) and returns 'succeeded' if movement is successful 
        
        Args:
            ent_name -- The name of the entity to move
            direction -- a three integer tuple with two zero entries and one +1 or -1,
                specifying a direction to attempt to move the entity in
        """
        
        # First we find the entity by name
        number = None
        ent = None
        for key,entity in self.entities.items():
            if entity.name == ent_name:
                ent = entity
                number = key
                break
        else:
            raise ValueError('Entity {} not found'.format(ent_name))
        
        i,j,k = ent.coords
        diri,dirj,dirk = direction
        newi,newj,newk = i+2*diri,j+2*dirj,k+2*dirk
        
        if self.maze[i+diri,j+dirj,k+dirk] == 0 and self.maze[newi,newj,newk] == 0:
            # If move successful, update entries in self.maze and ent
            self.maze[newi,newj,newk] = number
            self.maze[i,j,k] = 0
            ent.coords = (newi,newj,newk)
            return ent
        elif self.maze[i+diri,j+dirj,k+dirk] == 1: 
            # If we collided with a wall
            # return a wall entity showing where we collided
            return MazeEntity('wall', 
                              self.display_options['wall_symbol'], 
                              (i+diri,j+dirj,k+dirk))
        else: # If we collided with another entity
            return self.entities[self.maze[newi,newj,newk]]

direction_mapping = {'w' : (0,-1,0),
                     'a' : (0,0,-1),
                     's' : (0,1,0),
                     'd' : (0,0,1)}

rotate_mapping = {'rw' : CubeIsometry('j -i k'),
                  'ra' : CubeIsometry('-k j i'),
                  'rs' : CubeIsometry('-j i k'),
                  'rd' : CubeIsometry('k j -i')}

def play_game(n):
    maze_game = MazeNav3D(n,display_options='med')
    player = MazeEntity('player',('/\\','\\/'),(1,1,1))
    goal = MazeEntity('goal',('OO','OO'), (-2,-2,-2))
    maze_game.add_entity(player)
    maze_game.add_entity(goal)
    
    print("Input one command at a time to move around.\n\
Commands are 'w' 'a' 's' 'd' for movement and 'rw' 'ra' 'rs' 'rd'\
to rotate the maze\nInput q or exit to quit.")
    while(True):
        maze_game.display_all()
        inp = input()
        if inp == 'q' or inp == 'exit':
            print('Goodbye')
            break
        if inp in direction_mapping:
            prev_coords = player.coords
            res = maze_game.attempt_move('player',direction_mapping[inp])
            if res.name == 'goal':
                print('You made it to the goal!  Amazing!')
                break
            elif res.name == 'wall':
                print('Ouch!  You hit a wall.')
            elif res.name[:4] == 'void':
                print('You were destroyed by the void!')
                break
            else:
                noreturn = MazeEntity('void' + str(maze_game.next_entity - 2),
                                      ('XX','XX'),
                                      prev_coords)
                maze_game.add_entity(noreturn)
        elif inp in rotate_mapping:
            maze_game.apply_isometry(rotate_mapping[inp])
            
if __name__ == "__main__":
    play_again = True
    while(play_again):
        play_game(3)
        print('Play again? (y/n)')
        inp = input()
        while(inp != 'y' and inp != 'n'):
            inp = input()
        if inp == 'n':
            play_again = False