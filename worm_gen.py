# Usage on command line:
#   First parameter = Number of desired worms
#   2nd & 3rd = min & max of Km
#   4th & 5th = min & max of N
#   6th & 7th = min & max of GJ_scale
#   8th & 9th = min & max of cells_num
#   10th & 11th = min & max of G_k 
#   12th & 13th = min & max of G_na
#   14th & 15th = min & max of G_cl
#   16th & 17th = min & max of GJ_diff_m

NUM_TRAITS = 9

import random
import sys
import fileinput
import string

class Worm:
    def __init__(self, trait_bounds):
        self.Km = random.randint(trait_bounds[1], trait_bounds[2])
        self.N = random.randint(trait_bounds[3], trait_bounds[4])
        self.Gj_scale = random.randint(trait_bounds[5], trait_bounds[6])
        self.num_cells = random.randint(trait_bounds[7], trait_bounds[8])
        self.G_k = random.randint(trait_bounds[9], trait_bounds[10])
        self.G_na = random.randint(trait_bounds[11], trait_bounds[12])
        self.G_cl = random.randint(trait_bounds[13], trait_bounds[14])
        self.Gj_diff_m = random.randint(trait_bounds[15], trait_bounds[16])
        self.time_to_stable = 0
        self.length_of_stable = 0
        self.grad_stren = 0
        self.fitness = 0

    def __lt__(self, other):
        if(self.fitness < other.fitness):
            return self
        else:
            return other
    # def __eq__(self, other):
    #     if(self.fitness == other.fitness):
    #         return self
    #     else:
    #         return other
    # def __hash__(self):
    #     return 

def gen_worms(trait_file):
    trait_bounds = get_traits(trait_file)
    trait_bounds = [int(i) for i in trait_bounds]
    num_worms = trait_bounds[0]
    random.seed()
    i = 0
    population = set()
    while(i < num_worms):
        worm = Worm(trait_bounds)
        population.add(worm)
        i = i + 1
    #print_pop(population, num_worms)
    return population, trait_bounds

def get_traits(trait_file):
    i = 1
    trait_bounds = []
    while (i <= (NUM_TRAITS * 2)):
        line = trait_file.readline()
        if i % 2 == 0:
            seperate_traits = line.split()
            if((i > 2) and (int(seperate_traits[0]) > int(seperate_traits[1]))):
                print("\nIncorrectly formatted input file.\nPlease keep original"+
                    " documetnation and ensure all lower bounds are less than"+
                    " higher bounds.\n")
                exit(1)
            trait_bounds.extend(seperate_traits) 
        i = i + 1
    return trait_bounds



