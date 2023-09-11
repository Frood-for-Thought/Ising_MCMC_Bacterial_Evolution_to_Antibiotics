

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

# This algorithm determines the neighbouring spins next to the spin at the position selected and calculates
# the necessary values for the energy equation.  It's made so that positions at the edges and corners are
# not connected to the opposite edges like a toroid shape, and instead the grid is finite.
def Detailed_Bal_Eqn_Val(array, row, col, N, shift_list=None):
    A = []
    B = []
    C = []
    D = []
    if shift_list is None:
        shift_list = [{"Label": "A", "Include": 0, "Matrix": A},
                      {"Label": "B", "Include": 0, "Matrix": B},
                      {"Label": "C", "Include": 0, "Matrix": C},
                      {"Label": "D", "Include": 0, "Matrix": D}]

    # Determine if the spin is on the left or right corner, or neither
    if col == 0:
        shift_list[0]["Include"] = 1
    elif col == N - 1:
        shift_list[1]["Include"] = 1
    else:
        shift_list[0]["Include"] = 1
        shift_list[1]["Include"] = 1
    # Determine if the spin is on the top or bottom corner, or neither
    if row == 0:
        shift_list[2]["Include"] = 1
    elif row == N - 1:
        shift_list[3]["Include"] = 1
    else:
        shift_list[2]["Include"] = 1
        shift_list[3]["Include"] = 1

    # Sort the dictionary so the matrices that will be calculated are at the beginning of the list
    sl = sorted(shift_list, key=lambda x: x["Include"] > 0, reverse=True)

    # Need to create tuple in reverse so that every element is counted as the mutable list is shortened,
    # (this would not be the case if the list went forwards)
    # The neighbouring spins are recorded by shifting the entire array so that the neighbouring values are at the
    # position, (row, col), where the center spin is located
    for i, j in reversed(list(enumerate(sl))):
        # The "Include" variable is 1 to include and 0 to not include
        if j["Include"] > 0:
            # Calculate and record information for Matrices A, B, C, or D
            # Calculate Matrix Shifts
            if j["Label"] == "A":
                # A = moves all elements to the right[0 1], replaces circshift(spin, [0 1]) in MATLAB
                A = np.roll(array, -1)
                j["Matrix"] = A[row][col]
            elif j["Label"] == "B":
                # B = moves all elements to the left[0 - 1], replaces circshift(spin, [0 - 1]) in MATLAB
                B = np.roll(array, 1)
                j["Matrix"] = B[row][col]
            elif j["Label"] == "C":
                # C = Shifts down all the elements, replaces circshift(spin, [1 0]), in MATLAB # numpy axis are the
                # directions along the rows and columns, with axis=0 being the rows and axis=1 being the columns
                C = np.roll(array, -1, axis=0)
                j["Matrix"] = C[row][col]
            elif j["Label"] == "D":
                # D = Shift up moves all the elements, replaces circshift(spin, [-1 0]) in MATLAB
                D = np.roll(array, 1, axis=0)
                j["Matrix"] = D[row][col]
        else:
            # Start with removing the matrices that are not used, sl.remove(sl[i])
            sl.pop()

    # sum all neighbouring spins (NS)
    neighbours = sum(mtx["Matrix"] for mtx in sl)
    # sum NS^2
    neighbours_sqrd = sum(np.square(mtx["Matrix"]) for mtx in sl)
    # sum (1 - NS)
    OneMinusNeighbour = sum(np.subtract(1, mtx["Matrix"]) for mtx in sl)
    # sum (1 - NS^2)
    death_neighbour = sum(np.subtract(1, np.square(mtx["Matrix"])) for mtx in sl)
    # Check to see if all neighbours are present to use the Food Affinity part of Master Equation
    Counter_Food_Affinity = 0
    for i in range(len(sl)):
        Counter_Food_Affinity += np.prod(sl[0]["Matrix"])
    return neighbours, neighbours_sqrd, OneMinusNeighbour, death_neighbour, Counter_Food_Affinity