import heapq
import math

PROP_WORMS_FAIL = 0.50

def population_fitness(population):
    pop_size = len(population)
    pop_heap = []
    i = 0
    #populationC = population.copy()

    #heapq.heappush(pop_heap, (worm.fitness, worm))
    while i < pop_size:
        worm = population.pop()
        if (type(worm) != tuple):
            worm.fitness = calculate_fitness(worm)
            pop_heap.append((worm.fitness, worm))
        if(type(worm) == tuple):
            worm[1].fitness = calculate_fitness(worm[1])
            pop_heap.append(worm)
        i = i + 1
    heapq.heapify(pop_heap)
    #print(pop_heap)
    return pop_heap


def calculate_fitness(worm):
    #fitness function needs more complexity, but this is a start
    #GET RID OF 'i' DUMBO
    fitness = worm.grad_stren
    return fitness

def purge_bad_worms(population):
    pop_size = len(population)
    num_to_kill = math.trunc(PROP_WORMS_FAIL * pop_size)
    while (num_to_kill > 0):
        heapq.heappop(population)
        num_to_kill = num_to_kill - 1
    return population



