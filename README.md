# WormEvolution
Evolutionary Algorithm for simulating Planaria using BITSEY.

Files developed for this project: main.py, bitsey.py, test_wormAvgs.py, worm_breeding.py, worm_fitness.py, worm_gen.py, worm_input.txt.

Bitsey (a simulation developed by Alexis Pietak and Joel Grodstein) is given a set of parametrs and calculates the fitness of each worm using BITSEY (Cell matrix voltage simulator).

Our algorithm takes the fittest 50% of the population and doubles this most-fit set by "breeding" them in pairs.
Each pair of worms produces two children with a mixture of their parents' traits original traits as well as randomly selected "mutated" traits. 

The goal of this project is to offer researchers a tool with which they can test to see which parameters produce the highest
fitness levels among a population of worms.
