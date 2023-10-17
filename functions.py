import numpy as np


class Ising_Functions:
    def __init__(self, array, pos, N):
        self.array = array
        self.pos = pos
        self.N = N
        self.row, self.col = self.fnd_pos()

    def fnd_pos(self):
        """
        This replaces the MATLAB built-in function: [row, col] = ind2sub(array, rand_element_selected).
        Given the array and the random element number in the array, find the row and column for said element the elements
        are numbered from left to right, top to bottom.
        :param array: Insert the array.
        :param pos: Random array element number inserted.
        :return: The row, col position for the array element number.
        """
        col = self.pos - 1  # start at column 0
        for i, j in enumerate(self.array):
            col -= len(j)
            if col < 0:
                col += len(j)
                row = i
                return row, col

    def Fitness_Eqn_Val(self, shift_list=None):
        """
        This algorithm determines the neighbouring spins next to the spin at the position selected and calculates
        the necessary values for the energy equation.  It's made so that positions at the edges and corners are
        not connected to the opposite edges like a toroid shape, and instead the grid is finite.  The function makes
        a new matrix which shifts the positions of the neighbouring spins to be located where the central spin is.
        :param array: Insert the array.
        :param row: The row located for the selected element.
        :param col: The col located for the selected element.
        :param N: The N x N size of "array".
        :param shift_list: A dictionary list for the elements neighbouring the central selected element at (row, col).
        Each element dictionary contains a label for the position (up, down, left, right), an Include integer (0 or 1),
        and a Matrix to shift the spin towards the position of the selected element (row, col).
        :return:
        """
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
        if self.col == 0:
            shift_list[0]["Include"] = 1
        elif self.col == self.N - 1:
            shift_list[1]["Include"] = 1
        else:
            shift_list[0]["Include"] = 1
            shift_list[1]["Include"] = 1
        # Determine if the spin is on the top or bottom corner, or neither
        if self.row == 0:
            shift_list[2]["Include"] = 1
        elif self.row == self.N - 1:
            shift_list[3]["Include"] = 1
        else:
            shift_list[2]["Include"] = 1
            shift_list[3]["Include"] = 1

        # Sort the dictionary so the matrices that will be calculated are at the beginning of the list
        sl = sorted(shift_list, key=lambda x: x["Include"] > 0, reverse=True)

        # Need to create tuple in reverse so that every element is counted as the mutable list is shortened
        # for each iteration.
        for i, j in reversed(list(enumerate(sl))):
            # The "Include" variable is 1 to include and 0 to not include
            if j["Include"] > 0:
                # The neighbouring spins are recorded by shifting the entire array so that the neighbouring values
                # are at the position, (row, col), where the center spin is located.
                # Calculate and record information for Matrices A, B, C, or D
                # Calculate Matrix Shifts
                if j["Label"] == "A":
                    # A = moves all elements to the right[0 1], replaces circshift(spin, [0 1]) in MATLAB
                    A = np.roll(self.array, -1)
                    j["Matrix"] = A[self.row][self.col]
                elif j["Label"] == "B":
                    # B = moves all elements to the left[0 - 1], replaces circshift(spin, [0 - 1]) in MATLAB
                    B = np.roll(self.array, 1)
                    j["Matrix"] = B[self.row][self.col]
                elif j["Label"] == "C":
                    # C = Shifts down all the elements, replaces circshift(spin, [1 0]), in MATLAB # numpy axis are the
                    # directions along the rows and columns, with axis=0 being the rows and axis=1 being the columns
                    C = np.roll(self.array, -1, axis=0)
                    j["Matrix"] = C[self.row][self.col]
                elif j["Label"] == "D":
                    # D = Shift up moves all the elements, replaces circshift(spin, [-1 0]) in MATLAB
                    D = np.roll(self.array, 1, axis=0)
                    j["Matrix"] = D[self.row][self.col]
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

    @staticmethod
    def partition_gillespie(prob_list, rand, R=None, Cond_list=None):
        """
        Use partition function to find which transition occurs in a manner similar to the 2nd stage of the
        Gillespie Algorithm. This function basically gets a list of varying length with probabilities as the elements,
        then loops through to calculate the rates of a reaction occurring.
        It then uses the rates to determine which condition occurred and marks one as True while the rest are False.
        :param prob_list: Insert in a list of all the probability of the allowed reactions.
        :param rand: A previously generated random number that is inserted.
        :param R: A list constructed to contain all the normalized rates or reactions taken from prob_list.
        :param Cond_list: List of booleans, with default = False.  When rand is in the normalized range for a rate
        of reaction to occur, the condition for that reaction will be set to True.
        :return: Cond_list
        """
        p_sum = sum(prob_list)
        # Set R and Cond_list as empty arrays.
        if R is None and Cond_list is None:
            R = []
            Cond_list = []
        for i, p in enumerate(prob_list):
            if i < 1:
                R.append(p / p_sum)
                Cond_list.append(False)
            else:
                R.append(R[i - 1] + p / p_sum)
                Cond_list.append(False)
        for r in range(len(R)):
            if r == 0 and rand <= R[r]:
                Cond_list[r] = True
            elif R[r - 1] < rand <= R[r]:
                Cond_list[r] = True
        return Cond_list

    @staticmethod
    def allow_transition_state(allow_fit, fit_list, spin=0):
        """
        A function to determine which transition is allowed.  The fit_list is selected depending on the previous
        state the spin was in, and is used to determine which transition is allowed next.
        :param spin: return which spin the next state will be.
        :param allow_fit = sum([i["bool"] for i in fit_list])
        :param fit_list: a list of dictionaries which each contain a boolean to determine if the fitness (energy) is
        negative, the probability of switching to that state if it is not favorable (due to detailed balance),
        and the final spin state if that transition is allowed.
        :return:
        """
        if allow_fit > 1:
            prob_rand = np.random.rand()
            # Use partition function to select which transition occurs.
            prob_list = [j["prob"] for j in fit_list]
            Rnd_Select = Ising_Functions.partition_gillespie(prob_list, prob_rand)
            for k, l in enumerate(Rnd_Select):
                if l:
                    spin = fit_list[k]["spin"]
        # Only one accepted.
        elif allow_fit > 0:
            for j in fit_list:
                if j["bool"]:
                    spin = j["spin"]
        # Probability of transition still occurring due to detailed balance.
        else:
            prob_rand_g = np.random.rand()
            # The probability of both is selected
            # Use transition probabilities from detailed balance to find which transition occurs.
            allow_prob = sum((prob_rand_g < i["prob"] for i in fit_list))
            # Both probabilities from detailed balance are selected due to being above the random number selected.
            if allow_prob > 1:
                # Reroll number again to randomize which transition occurs.
                prob_rand = np.random.rand()
                # Use partition function to find which transition occurs.
                prob_list = [j["prob"] for j in fit_list]
                Rnd_Select = Ising_Functions.partition_gillespie(prob_list, prob_rand)
                for k, l in enumerate(Rnd_Select):
                    if l:
                        spin = fit_list[k]["spin"]
            # Probability of 0 --> 1 selected or Probability of 0 --> -1 selected
            else:
                for j in fit_list:
                    if prob_rand_g < j["prob"]:
                        spin = j["spin"]
        return spin
