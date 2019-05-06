import heapq
import math

PROP_WORMS_FAIL = 0.50

def population_fitness(population):
    pop_size = len(population)
    pop_heap = []
    i = 0
    #heapq.heappush(pop_heap, (worm.fitness, worm))
    while i < pop_size:
        worm = population.pop()
        worm.fitness = calculate_fitness(worm, i)
        pop_heap.append((worm.fitness, worm))
        i = i + 1
    heapq.heapify(pop_heap)
    #print(pop_heap)
    return pop_heap


def calculate_fitness(worm, i):
    #fitness function needs more complexity, but this is a start
    #GET RID OF 'i' DUMBO
    fitness = worm.grad_stren + worm.length_of_stable - worm.time_to_stable + i
    return fitness

def purge_bad_worms(population):
    pop_size = len(population)
    num_to_kill = math.trunc(PROP_WORMS_FAIL * pop_size)
    while (num_to_kill > 0):
        heapq.heappop(population)
        num_to_kill = num_to_kill - 1
    return population



