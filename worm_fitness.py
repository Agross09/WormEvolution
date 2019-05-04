import heapq


def population_fitness(population):
    pop_size = len(population)
    pop_heap = []
    i = 0
    #heapq.heappush(pop_heap, (worm.fitness, worm))
    while i < pop_size:
        worm = population.pop()
        worm.fitness = calculate_fitness(worm, i)
        heapq.heappush(pop_heap, (worm.fitness, worm))
        i = i + 1
    return pop_heap




def calculate_fitness(worm, i):
    fitness = worm.grad_stren + worm.length_of_stable - worm.time_to_stable
    return fitness
