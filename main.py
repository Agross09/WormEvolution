import worm_gen 
import worm_fitness as fit
import worm_breeding as breed
from heapq import merge, heapify, heappop
import bitsey as bit
import sys
import test_wormAvgs as wormAvg

# Notes: give worm boolean to see if tested with bitsey. DO NOT want to run worm through bitsey more than once.

#Determines number of generations to be created
NUM_GENS = 8

def main():
    trait_file = open(sys.argv[1], "r")
    population, trait_bounds = worm_gen.gen_worms(trait_file)

    #Simulation time for each worm
    time = int(input("How long do you want the worms to run?"))

    i = 0
    population = bitsey(population,time)
    #main loop. Each iteration is one generation
    while i < NUM_GENS:
        print("Gen {} -----------------".format(i))

        #Runs each worm in population through bitsey to determine gradient 
        #    strengths
        
        #print(population)
        
        # if i == 1:
        #     #prints gradient strength of original worms
        #     print("Parents Grads: ")
        #     wormResults(population)
        #     print("--------------")

        #determines fitness for each worm in population
        sorted_pop = fit.population_fitness(population)

        #removes least fit 50% of worms
        best_worms_pop = fit.purge_bad_worms(sorted_pop)

        #breeds worms and returns a list of children
        child_list = breed.breed_worms(best_worms_pop, trait_bounds)
        sorted_kids = bitsey(child_list, time)
        # an unordered list of new worms without fitness levels - not BITSEYed yet.
                           # Doubles population
        # Make heap of child worms with fitness levels gotten through BITSEY:
        child_pop = fit.population_fitness(sorted_kids)
        # Pass child_pop and the best_worms_pop into a function that merges them into a single heap:
        population = mergeHeaps(best_worms_pop, child_pop)


        i = i + 1
        
        wormResults(population)
        print("-------------------------\n")
        
    #prints average KM, N, Gj_scale values. Not important, just from old tests
    ##wormAvg.average_Km(population)
    #wormAvg.average_N(population)
    #wormAvg.average_Gj_scale(population)

    #prints gradient strength of each worm
    
    trait_file.close()
    return

def wormResults(pop):
    population = pop.copy()
    pop_size = len(population)
    to_print = list()
    i = 0
    tot = 0
    while i < pop_size:
        worm = heappop(population)
        to_print.append(worm)
        i = i + 1
    i = pop_size
    while i > 0:
        worm = to_print[i - 1]
        if (type(worm) == tuple):
            #print("Grad Stren: {}".format(worm[1].fitness))
            tot = tot + worm[1].fitness
        else:
            print("nGrad Stren: {}".format(worm.fitness))
        i = i - 1
        avg = tot / pop_size
    print("Average Fitness: {}".format(avg))


def bitsey(population,time):
    pop_size = len(population)
    #popC = population.copy()
    i = 0
    #popC = population.copy()
    if type(population) == set:
        #print("Set")
        popC = set()
        while i < pop_size:
            worm_to_work = population.pop()
            worm_to_work = bit.setup_and_sim(worm_to_work, time)
            popC.add(worm_to_work)
            i = i + 1

    if type(population) == list:
       # print("list")
        popC = list()
        while i < pop_size:
            worm_to_work = population.pop(0)
            worm_to_work.fitness = bit.setup_and_sim(worm_to_work, time).fitness
            popC.append(worm_to_work)
            i = i + 1
    return popC


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