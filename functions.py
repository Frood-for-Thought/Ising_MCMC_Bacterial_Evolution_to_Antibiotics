# This replaces the MATLAB built-in function: [row, col] = ind2sub(array, rand_element_selected)
# Given the array and the random element number in the array, find the row and column for said element
# the elements are numbered from left to right, top to bottom
def fnd_pos(array, pos):
    col = pos - 1  # start at column 0
    for i, j in enumerate(array):
        col -= len(j)
        if col < 0:
            col += len(j)
            row = i
            return row, col

# Construct initial Gillespie rate probability function

# Construct function to calculate equation values to be inserted at each position selected.