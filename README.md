# Ising_MCMC_Bacterial_Evolution_to_Antibiotics
A Markov Chain Monte Carlo Ising Model for Bacterial Evolution. Using the mathematical model in ferromagnetism in statistical mechanics called the Ising Model as a basis 
in order to create a new program to simulate how an antibiotic gradient affects the rate of mutation of a bacterial population.
The program uses a Hamiltonian equation to govern which direction a particular microstate, (a bacterial population in a deme), will progress towards 
in order to describe bacterial populations migrating along food and drug gradients.  There are several factors descibed in the Hamiltonian which affect the evolutionary 
rate, the growth rate, death rate, and rate of migration. 


Unlike the classic Ising model, the spin values (or demes) for this system can have three different values:
- A deme which is fixed with a population of wild type bacteria is given a +1. 
- A deme which is absent of any bacteria is given a value of 0. 
- A deme which represents a system which has been fixed by mutant bacteria is given a value of -1.


The Hamiltonian can be described as:

```math
E = -J/2 \sum_{i \neq j}^N S_i S_j + J_d \sum_{i=1}^L (S_i)^2 - \sum_{i \neq j}^N J_f/2*(S_i)^2 (S_j)^2 + \sum_{i=1}^L J_c*S_i
```
```math
- A*[\sum_{i \neq j}^N exp(J_fmax/kT + J_f/kT + 0.095)(1-S_i)(1-S_j)][\prod_{i = 1}^L (S_i)^2]
```
- N are nearest neighbours, (N = 4)
- L are all neighbouring and central spins, (L = 5)
- J is the exchange constant to support similar bacterial populations.
- ```math
J_d
``` is the bacterial death term to shift populated areas to empty (+/- 1 to 0).
- J_f is the food term that's a variable dependent on the position and supports increasing popoulation so long as the food supply is sufficient.
- J_c is the antibiotic term and the energy changes in relation to the bacterial present.  
The term supports the growth of mutant bacteria while discouraging the growth of wild-type bacteria.
- A is food affinity term to support mutant growth by giving it a greater affinity than wild-type.
- kT is the temperature constant in thermodynamics, but in this case represents the stochastic drive of a gene, (a larger kT causing greater genetic variation).

## Citation
The work is made using the thesis "Simulated accelerated evolution by modeling the rapid fixation of bacteria along an antibiotic gradient":
https://knowledgecommons.lakeheadu.ca/handle/2453/4660
