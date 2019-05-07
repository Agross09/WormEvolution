import worm_gen 
import worm_fitness as fit
import worm_breeding as breed
from heapq import merge, heapify
#import main as bitsey
import sys

# Notes: give worm boolean to see if tested with bitsey. DO NOT want to run worm through bitsey more than once.

NUM_GENS = 10000

def main():
    trait_file = open(sys.argv[1], "r")
    population, trait_bounds = worm_gen.gen_worms(trait_file)
    i = 0
    print("hey")
    while i < NUM_GENS:
        if i == 1:
            print("Parents Avgs: ")
            parents_avg(population)
            print("--------------")
        sorted_pop = fit.population_fitness(population)
        best_worms_pop = fit.purge_bad_worms(sorted_pop)
        child_list = breed.breed_worms(best_worms_pop, trait_bounds)
        # an unordered list of new worms without fitness levels - not BITSEYed yet.
                           # Doubles population
        # Make heap of child worms with fitness levels gotten through BITSEY:
        child_pop = fit.population_fitness(child_list)
        # Pass child_pop and the best_worms_pop into a function that merges them into a single heap:
        population = mergeHeaps(best_worms_pop, child_pop)
        i = i + 1
    average_Km(population)
    average_N(population)
    average_Gj_scale(population)
    trait_file.close()
    return

def parents_avg(population):
    average_Km(population)
    average_N(population)
    average_Gj_scale(population)

def average_Km(end_of_round_pop):
    worm_number = len(end_of_round_pop)
    denominator = len(end_of_round_pop)
    total = 0
    while worm_number > 0:
        total = end_of_round_pop[worm_number - 1][1].Km + total
        worm_number = worm_number - 1
    avg = total / denominator
    print(" Avg KM: {}".format(avg))


def average_N(end_of_round_pop):
    worm_number = len(end_of_round_pop)
    denominator = len(end_of_round_pop)
    total = 0
    while worm_number > 0:
        total = end_of_round_pop[worm_number - 1][1].N + total
        worm_number = worm_number - 1
    avg = total / denominator
    print(" Avg N: {}".format(avg))

def average_Gj_scale(end_of_round_pop):
    worm_number = len(end_of_round_pop)
    denominator = len(end_of_round_pop)
    total = 0
    while worm_number > 0:
        total = end_of_round_pop[worm_number - 1][1].Gj_scale + total
        worm_number = worm_number - 1
    avg = total / denominator
    print(" Avg GJ: {}".format(avg))


def print_heap(population, num_worms):
    i = 0
    print("\n")
    while(i < num_worms):
        to_print = population[i]
        print("Worm Fitness = " + str(to_print[1].fitness))
        print("Km = " + str(to_print[1].Km))
        print("N = " + str(to_print[1].N))
        print("Gj_scale = " + str(to_print[1].Gj_scale))
        print("num_cells = " + str(to_print[1].num_cells))
        print("G_k = " + str(to_print[1].G_k))
        print("G_k = " + str(to_print[1].G_na))
        print("G_cl = " + str(to_print[1].G_cl))
        print("Gj_diff_m = " + str(to_print[1].Gj_diff_m))
        print("value")
        print("\n")
        i = i + 1

def mergeHeaps(heap1, heap2):
    return list(merge(heap1, heap2))
 

main()