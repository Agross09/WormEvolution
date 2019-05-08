import worm_gen 
import worm_fitness as fit
import worm_breeding as breed
from heapq import merge, heapify, heappop
import bitsey as bit
import sys
import test_wormAvgs as wormAvg

# Notes: give worm boolean to see if tested with bitsey. DO NOT want to run worm through bitsey more than once.

NUM_GENS = 2

def main():
    trait_file = open(sys.argv[1], "r")
    population, trait_bounds = worm_gen.gen_worms(trait_file)
    i = 0
    time = int(input("How long do you want the worms to run?"))
    population_list = changeSetToListofTuples(population)
    population_list = bitsey(population_list,time)
    while i < NUM_GENS:
        print("Gen {}".format(i))
        if i == 1:
            print("Parents Grads: ")
            wormAvg.parents_avg(population_list)
            wormResults(population_list)
            print("--------------")
        sorted_pop = fit.population_fitness(population_list)
        best_worms_pop = fit.purge_bad_worms(sorted_pop)
        child_list = breed.breed_worms(best_worms_pop, trait_bounds)
        # an unordered list of new worms without fitness levels - not BITSEYed yet.
                           # Doubles population
        # Make heap of child worms with fitness levels gotten through BITSEY:
        child_list = bitsey(child_list,time)
        child_pop = fit.population_fitness(child_list)
        # Pass child_pop and the best_worms_pop into a function that merges them into a single heap:
        population_list = mergeHeaps(best_worms_pop, child_pop)


        i = i + 1
    wormAvg.average_Km(population_list)
    wormAvg.average_N(population_list)
    wormAvg.average_Gj_scale(population_list)
   # heapify(population)
    wormResults(population_list)
    trait_file.close()
    return

def wormResults(pop):
    population = pop.copy()
    print(population)
    pop_size = len(population)
    to_print = list()
    i = 0
    while i < pop_size:
        worm = heappop(population)
        to_print.append(worm)
        i = i + 1
    i = pop_size
    while i > 0:
        worm = to_print[i - 1]
        print("Grad Stren: {}".format(worm[0]))
        i = i - 1


def bitsey(population,time):
    pop_size = len(population)
    i = 0
    popC = population.copy()
    if type(popC) == list:
        while i < pop_size:
            worm_to_work = popC[i]
            print("BEFORE Bitsey: {}".format(worm_to_work[1].grad_stren))
            bit.setup_and_sim(worm_to_work[1], time)
            print("AFTER Bitsey: {}".format(worm_to_work[1].grad_stren))
            popC[i] = worm_to_work
            i = i + 1
        return popC

def changeSetToListofTuples(population):
    if type(population) == set:
        i = 0
        size_of_set = len(population)
        pop_list = list()
        fitness = 0
        while i < size_of_set:
            worm_tup = (fitness, population.pop())
            pop_list.append(worm_tup)
            i = i + 1
        return pop_list

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