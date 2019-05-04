import worm_gen 
import worm_fitness as fit
#import main as bitsey
import sys



def main():
    trait_file = open(sys.argv[1], "r")
    population = worm_gen.gen_worms(trait_file)
    sorted_pop = fit.population_fitness(population)
    print_heap(sorted_pop, len(sorted_pop))
    trait_file.close()
    return



def print_heap(population, num_worms):
    i = 0
    print("\n")
    while(i < num_worms):
        to_print = population[i]
        print("Worm Fitness = " + str(to_print[1].fitness))
        print("Km = " + str(to_print[1].Km))
        print("N = " + str(to_print[1].N))
        print("GJ_scale = " + str(to_print[1].GJ_scale))
        print("num_cells = " + str(to_print[1].num_cells))
        print("G_k = " + str(to_print[1].G_k))
        print("G_k = " + str(to_print[1].G_na))
        print("G_cl = " + str(to_print[1].G_cl))
        print("GJ_diff_m = " + str(to_print[1].GJ_diff_m))
        print("value")
        print("\n")
        i = i + 1


main()