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
        self.GJ_scale = random.randint(trait_bounds[5], trait_bounds[6])
        self.num_cells = random.randint(trait_bounds[7], trait_bounds[8])
        self.G_k = random.randint(trait_bounds[9], trait_bounds[10])
        self.G_na = random.randint(trait_bounds[11], trait_bounds[12])
        self.G_cl = random.randint(trait_bounds[13], trait_bounds[14])
        self.GJ_diff_m = random.randint(trait_bounds[15], trait_bounds[16])

        #working set of parameters
        # self.Km = random.randint(1,1)
        # self.N = random.randint(10,10)
        # self.GJ_scale = random.randint(1, 1)
        # self.num_cells = random.randint(5,5)
        # self.G_k = random.randint(trait_bounds[9], trait_bounds[10])
        # self.G_na = random.randint(trait_bounds[11], trait_bounds[12])
        # self.G_cl = random.randint(trait_bounds[13], trait_bounds[14])
        # self.GJ_diff_m = 1e-14

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
    
    return population

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



def print_pop(population, num_worms):
    i = 0
    print("\n")
    print_pop = population.copy()
    while(i < num_worms):
        to_print = print_pop.pop()
        print("Worm #" + str(i))
        print("Km = " + str(to_print.Km))
        print("N = " + str(to_print.N))
        print("GJ_scale = " + str(to_print.GJ_scale))
        print("num_cells = " + str(to_print.num_cells))
        print("G_k = " + str(to_print.G_k))
        print("G_k = " + str(to_print.G_na))
        print("G_cl = " + str(to_print.G_cl))
        print("GJ_diff_m = " + str(to_print.GJ_diff_m))
        print("\n")
        i = i + 1
worm_gen.py
Displaying worm_gen.py.