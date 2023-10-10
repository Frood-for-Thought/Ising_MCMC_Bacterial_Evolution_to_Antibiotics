import random
import numpy as np
import matplotlib.pyplot as plt

# A MARKOV CHAIN MONTE CARLO (MCMC) METHOD TO DESCRIBE BACTERIAL EVOLUTION USING THE ISING MODEL

# All the parameters are set up for N = 100
kT = 2
N = 100

# Check to make sure "numSpinsPerDim" is even.

Even = divmod(N, 2)
if Even[1] != 0:  # Odd number
    print('Cannot use this value. N has to be even.')
    exit()

# Set up a zero array to contain random spin grid
spin = [[0 for col in range(N)] for row in range(N)]
# 4 spins in the center given value "1"
HalfPos = int(N / 2) - 1
spin[HalfPos][HalfPos] = 1
spin[HalfPos + 1][HalfPos] = 1
spin[HalfPos][HalfPos + 1] = 1
spin[HalfPos + 1][HalfPos + 1] = 1

spin[0][0] = 1
spin[0][1] = 1
spin[1][0] = 1
spin[1][1] = 1

# Reformat the array so that it is numpy friendly
spin = np.array(spin)

ColMax = int((N / 2) - 1)
RowMax = int((N / 2) - 1)
ColMin = 0
RowMin = 0

# FOOD GRID
# Set up an array to contain the food gradient
Food_Function = [[0 for col in range(N)] for row in range(N)]
for Food_ElRow in range(N):
    for Food_ElCol in range(N):
        if Food_ElCol <= ColMax and Food_ElRow <= RowMax:
            Food_Function[Food_ElRow][Food_ElCol] = max((ColMax + 1 - Food_ElCol), (RowMax + 1 - Food_ElRow))
        elif Food_ElCol > ColMax and RowMax >= Food_ElRow:
            Food_Function[Food_ElRow][Food_ElCol] = max((Food_ElCol - ColMax), (RowMax + 1 - Food_ElRow))
        elif Food_ElCol <= ColMax and RowMax < Food_ElRow:
            Food_Function[Food_ElRow][Food_ElCol] = max((ColMax + 1 - Food_ElCol), (Food_ElRow - RowMax))
        elif Food_ElCol > ColMax and Food_ElRow > RowMax:
            Food_Function[Food_ElRow][Food_ElCol] = max((Food_ElCol - ColMax), (Food_ElRow - RowMax))
Food_Function = np.array(Food_Function)
# Further adjustments to the food gradient.
Food_Func = [[0 for col in range(N)] for row in range(N)]
for Food_ElRow in range(N):
    Food_Func[Food_ElRow] = list(map(lambda i: 3 + i*(2/(N/2)), Food_Function[Food_ElRow]))
Food_Func = np.array(Food_Func)

# CIPRO GRID (ANTIBIOTIC)
# Set up an array to contain the cipro gradient
Cipro_Function = [[0 for col in range(N)] for row in range(N)]
El_Col_Sum = 0
for Cipro_ElRow in range(N):
    for Cipro_ElCol in range(N):
        El_Col_Sum = Cipro_ElRow + Cipro_ElCol
        # Top Right Square
        if El_Col_Sum >= (N - 1) and Cipro_ElCol > N / 2 >= Cipro_ElRow:
            Cipro_Function[Cipro_ElRow][Cipro_ElCol] = max((Cipro_ElCol - ColMax), (Cipro_ElRow - RowMax))
        # Lower Left Square
        elif El_Col_Sum >= (N - 1) and Cipro_ElCol <= N / 2 < Cipro_ElRow:
            Cipro_Function[Cipro_ElRow][Cipro_ElCol] = max((Cipro_ElCol - ColMax), (Cipro_ElRow - RowMax))
        # Lower Right Square
        elif El_Col_Sum >= (N - 1) and Cipro_ElCol > N / 2 < Cipro_ElRow:
            Cipro_Function[Cipro_ElRow][Cipro_ElCol] = max((Cipro_ElCol - ColMax), (Cipro_ElRow - RowMax))
Cipro_Function = np.array(Cipro_Function)

# Further adjustments to the antibiotic gradient.
Cipro_Func = [[0 for col in range(N)] for row in range(N)]
for Cipro_ElRow in range(N):
    Cipro_Func[Cipro_ElRow] = list(map(lambda i: i * (2.115/(N/2)), Cipro_Function[Cipro_ElRow]))
Cipro_Func = np.array(Cipro_Func)
print(Cipro_Func)

# # Show a figure of the Food and Cipro gradients
# fig, axis = plt.subplots(1, 2)
# im1 = axis[0].imshow(Food_Func)
# axis[0].set_title("Food Gradient")
# fig.colorbar(im1, ax=axis[0])
# im2 = axis[1].imshow(Cipro_Func)
# fig.colorbar(im2, ax=axis[1])
# axis[1].set_title("Antibiotic Gradient")
# plt.show()

#   INITIALIZE THE ITERATION LOOP
# np.size(spin) counts the number of elements in the array "spin"
spin_array_size = np.size(spin)
print(f"The spin array has {spin_array_size} elements\n")
numIters = 200 * spin_array_size

# for curr_iter in range(numIters):
for curr_iter in range(5):

    # Pick a random spin in the array "spin"
    linearIndex = np.random.randint(1, spin_array_size)
    # print(f"Random Element Selected Number = {linearIndex}")

    # Use the function to find the position of the random spin
    spin_pos = fnd_pos(spin, linearIndex)
    row = spin_pos[0]
    col = spin_pos[1]

    # If there are no spins then the program moves on to the next iteration
    if spin[row][col] == 0:
        continue

    # GET VALUES FOR THE FITNESS EQUATION
    get_values = Fitness_Eqn_Val(spin, row, col, N)
    neighbours = get_values[0]
    # print(f"neighbours = {neighbours}")
    neighbours_sqrd = get_values[1]
    # print(f"neighbours_sqrd = {neighbours_sqrd}")
    one_minus_neighbour = get_values[2]
    # print(f"OneMinusNeighbour = {one_minus_neighbour}")
    neighbour_spin_product = get_values[4]
    # print(f"The product of the spin and its neighbours = {get_values[4]}")

    # If there are no spins then the program moves on to the next iteration
    if spin[row][col] == neighbours == 0:
        continue

    # EXCHANGE INTERACTION VALUES
    J = 10
    Jd = 40
    Jf = Food_Function[row][col]
    Jc = Cipro_Func[row][col]
    Jf_max = Food_Function[0][0]
    