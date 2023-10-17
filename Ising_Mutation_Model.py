import random
import numpy as np
import matplotlib.pyplot as plt

# A MARKOV CHAIN MONTE CARLO (MCMC) METHOD TO DESCRIBE BACTERIAL EVOLUTION USING THE ISING MODEL

# All the parameters are set up for N = 100.
kT = 0.2
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

# Import the Ising Model Functions Class
from functions import Ising_Functions

#   INITIALIZE THE ITERATION LOOP
# np.size(spin) counts the number of elements in the array "spin"
spin_array_size = np.size(spin)
print(f"The spin array has {spin_array_size} elements\n")
numIters = 200 * spin_array_size

for curr_iter in range(numIters):

    # Pick a random spin in the array "spin"
    linearIndex = np.random.randint(1, spin_array_size)

    # Use the function to find the position of the random spin
    row, col = Ising_Functions(spin, linearIndex, N).fnd_pos()

    # GET VALUES FOR THE FITNESS EQUATION
    get_values = Ising_Functions(spin, linearIndex, N).Fitness_Eqn_Val()
    neighbours = get_values[0]
    neighbours_sqrd = get_values[1]
    one_minus_neighbour = get_values[2]
    neighbour_spin_product = get_values[4]

    # If there are no spins then the program moves on to the next iteration
    if spin[row][col] == 0 and neighbours == 0:
        continue

    # EXCHANGE INTERACTION VALUES
    J = 1.86
    Jd = 5.95
    Jf = Food_Func[row][col]
    Jc = Cipro_Func[row][col]
    Jf_max = Food_Function[0][0]

    # SPIN SWITCH

    if spin[row][col] != 0:  # The selected spin is either -1 or 1

        # Death: spin = 1/-1 --> 0
        dE_1m1_0 = (J * spin[row][col] * neighbours) - Jd + (Jf * neighbours_sqrd) - (Jc * spin[row][col])
        - (J / 2) * (1 - spin[row][col]) * one_minus_neighbour * neighbour_spin_product * np.exp(-Jf_max + Jf + 0.095)

        # Spin Flip: spin = 1 --> -1, or spin = -1 --> 1
        dE_flip = (2 * J * spin[row][col] * neighbours) - (2 * Jc * spin[row][col])  # Putting 2*(J/2) below for clarity
        + 2 * (J / 2) * spin[row][col] * one_minus_neighbour * neighbour_spin_product * np.exp(-Jf_max + Jf + 0.095)

        # The probability of dying
        prob_d = np.exp(-dE_1m1_0 / kT)
        # The probability of flipping
        prob_flip = np.exp(-dE_flip / kT)

        # The fit_list is a list of dictionaries which each contain a boolean to determine if the fitness (energy) is
        # negative, the probability of switching to that state if it is not favorable (due to detailed balance),
        # and the final spin state if that transition is allowed.
        fit_list = [
            {"bool": dE_1m1_0 < 0, "prob": prob_d, "spin": 0},
            {"bool": dE_flip < 0, "prob": prob_flip, "spin": -spin[row][col]}]
        allow_fit = sum([i["bool"] for i in fit_list])

        # Select which next state the current spin is allowed to transition towards.
        spin[row][col] = Ising_Functions.allow_transition_state(allow_fit, fit_list, spin[row][col])

    else:  # The selected spin is 0
        # Growth: spin = 0 --> 1
        dE_0_1 = -(J * neighbours) + Jd - (Jf * neighbours_sqrd) + Jc
        + (J / 2) * one_minus_neighbour * neighbour_spin_product * np.exp(-Jf_max + Jf + 0.095)

        # Growth: spin = 0 --> -1
        dE_0_m1 = (J * neighbours) + Jd - (Jf * neighbours_sqrd) - Jc
        - (J / 2) * one_minus_neighbour * neighbour_spin_product * np.exp(-Jf_max + Jf + 0.095)

        # The probability of growing to 1
        prob_g_1 = np.exp(-(dE_0_1 / kT))
        # The probability of growing to -1
        prob_g_m1 = np.exp(-(dE_0_m1 / kT))

        # The fit_list is a list of dictionaries which each contain a boolean to determine if the fitness (energy) is
        # negative, the probability of switching to that state if it is not favorable (due to detailed balance),
        # and the final spin state if that transition is allowed.
        fit_list = [
            {"bool": dE_0_1 < 0, "prob": prob_g_1, "spin": 1},
            {"bool": dE_0_m1 < 0, "prob": prob_g_m1, "spin": -1}]
        allow_fit = sum([i["bool"] for i in fit_list])

        # Select which next state the current spin is allowed to transition towards.
        spin[row][col] = Ising_Functions.allow_transition_state(allow_fit, fit_list, spin[row][col])

    Display_Figure = divmod(curr_iter, 100)
    if Display_Figure[1] == 0:
        plt.clf()
        plt.imshow(spin)
        plt.colorbar(mappable=None, cax=None)
        plt.pause(0.05)
        plt.show(block=False)

    if curr_iter == numIters:
        print(spin)

# fig, axis = plt.subplots(1, 2)
# im1 = axis[0].imshow(Food_Func)
# axis[0].set_title("Food Gradient")
# fig.colorbar(im1, ax=axis[0])
# im2 = axis[1].imshow(Cipro_Func)
# fig.colorbar(im2, ax=axis[1])
# axis[1].set_title("Antibiotic Gradient")
# plt.show()
