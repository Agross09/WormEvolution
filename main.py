import worm_gen 
import worm_fitness as fit
import worm_breeding as breed
from heapq import merge, heapify
import bitsey as bit
import sys
import test_wormAvgs as wormAvg

# Notes: give worm boolean to see if tested with bitsey. DO NOT want to run worm through bitsey more than once.

NUM_GENS = 10

def main():
    trait_file = open(sys.argv[1], "r")
    population, trait_bounds = worm_gen.gen_worms(trait_file)
    i = 0
    while i < NUM_GENS:
        print("Gen {}".format(i))
        bitsey(population)
        if i == 1:
            print("Parents Avgs: ")
            wormAvg.parents_avg(population)
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
    wormAvg.average_Km(population)
    wormAvg.average_N(population)
    wormAvg.average_Gj_scale(population)
    trait_file.close()
    return

def bitsey(population):
    pop_size = len(population)
    i = 0
    popC = population.copy()
    if type(population) == set:
        while i < pop_size:
            worm_to_work = popC.pop()
            bit.setup_and_sim(worm_to_work, 10)
            i = i + 1
    if type(population) == list:
        while i < pop_size:
            worm_to_work = popC[i]
            bit.setup_and_sim(worm_to_work[1], 10)
            i = i + 1


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