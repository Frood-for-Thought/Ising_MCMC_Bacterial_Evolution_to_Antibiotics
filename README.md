# Ising_MCMC_Bacterial_Evolution_to_Antibiotics
A Markov Chain Monte Carlo Ising Model for Bacterial Evolution. Using the mathematical model in ferromagnetism in statistical mechanics called the Ising Model as a basis 
in order to create a new program to simulate how an antibiotic gradient affects the rate of mutation of a bacterial population.
The program uses a Hamiltonian equation to govern which direction a particular microstate, (a bacterial population in a deme), will progress towards 
in order to describe bacterial populations migrating along food and drug gradients.  There are several factors descibed in the Hamiltonian which affect the evolutionary 
rate, the growth rate, death rate, and rate of migration.

```math
E = -J/2 \sum_{i!=j}^N S_i S_j + J_d \sum_{i=1}^A (S_i)^2 - \sum_{i!=j}^N J_f/2*(S_i)^2 (S_j)^2 + \sum_{i=1}^A J_c*S_i
```
- N are nearest neighbours, (N = 4)
- A are all neighbouring and central spins, (A = 5)

Unlike the classic Ising model, the spin values (or demes) for this system can have three different values:
- A deme which is fixed with a population of wild type bacteria is given a +1. 
- A deme which is absent of any bacteria is given a value of 0. 
- A deme which represents a system which has been fixed by mutant bacteria is given a value of -1.

## Citation
The work is made using the thesis "Simulated accelerated evolution by modeling the rapid fixation of bacteria along an antibiotic gradient":
https://knowledgecommons.lakeheadu.ca/handle/2453/4660
