# Maze3D
A simple ascii game in which the player navigates a 3 dimensional maze to search for the goal

Now includes GUI interface with slightly simplified game.

TODO: 

1) Fix duplicate name checking in MazeNav3D.add_entity().  Possibly recode to allow duplicate names.  Maybe looking up objects by name is not necessary.  Better to have the add_entity function just return the code that can be used to look up the object.
2) Mess with main maze logic to make it more friendly to interfacing with GUI
3) Make GUI interface look nicer, and figure out how to get spinny arrows as the rotate buttons instead of text.
4) Possibly highlight the squares that will be moved to the top layer after rotating when the user hovers over a rotate button.
5) Add option for "non blocking" entity that can be passed through.  Maybe just store these in a list.  Then they can be moved around separately from the main maze objects.
