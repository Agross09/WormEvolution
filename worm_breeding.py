import worm_gen as worm
import heapq
import random

#Number of Kids per worm
NUM_KIDS = 2

#Chance to mutate a trait (1 / denom)
DENOM_CHANCE_TO_MUTATE = 2

#Worm Attribute Constants Do not affect data, just for readability)
KM = 1
N = 2
GJ_SCALE = 3
NUM_CELLS = 4
G_K = 5
G_NA = 6
G_CL = 7
GJ_DIFF_M = 8
NUM_TRAITS = 8

def breed_worms(population, trait_bounds):
    num_worms = len(population)
    i = num_worms/2 
    num_kids = NUM_KIDS
    worm_children_pop = list()
    population_copy = population.copy()
    while i > 0:
        worm1_with_fitness = heapq.heappop(population_copy)
        worm2_with_fitness = heapq.heappop(population_copy)
        worm1 = worm1_with_fitness[1]
        worm2 = worm2_with_fitness[1]
        while num_kids > 0:
            worm_baby = worm.Worm(trait_bounds)
            num_kids = num_kids - 1
            worm_baby = make_baby(worm1, worm2, worm_baby)
            worm_children_pop.append(worm_baby)
        num_kids = NUM_KIDS  
        i = i - 1
    return worm_children_pop


def make_baby(worm1, worm2, worm_baby):
    i = NUM_TRAITS
    from_worm1 = True
    
    
    # For each trait, we randomly choose from which parent we get the attribute.
    while i > 0:
        from_worm1 = random.getrandbits(1)
        chance_mutated = random.randint(1,DENOM_CHANCE_TO_MUTATE) # 1 in some number chance of mutation
        if chance_mutated != 1:
            if i == KM:
                if from_worm1 == True:
                    worm_baby.Km = worm1.Km
                else:
                    worm_baby.Km = worm2.Km
            if i == N:
                if from_worm1 == True:
                    worm_baby.N = worm1.N
                else:
                    worm_baby.N = worm2.N
            if i == GJ_SCALE:
                if from_worm1 == True:
                    worm_baby.Gj_scale = worm1.Gj_scale
                else:
                    worm_baby.Gj_scale = worm2.Gj_scale
            if i == NUM_CELLS:
                if from_worm1 == True:
                    worm_baby.num_cells = worm1.num_cells
                else:
                    worm_baby.num_cells = worm2.num_cells
            if i == G_K:
                if from_worm1 == True:
                    worm_baby.G_k = worm1.G_k
                else:
                    worm_baby.G_k = worm2.G_k
            if i == G_NA:
                if from_worm1 == True:
                    worm_baby.G_na = worm1.G_na
                else:
                    worm_baby.G_na = worm2.G_na
            if i == G_CL:
                if from_worm1 == True:
                    worm_baby.G_cl = worm1.G_cl
                else:
                    worm_baby.G_cl = worm2.G_cl
            if i == GJ_DIFF_M:
                if from_worm1 == True:
                    worm_baby.Gj_diff_m = worm1.Gj_diff_m
                else:
                    worm_baby.Gj_diff_m = worm2.Gj_diff_m

            # Iterate down through the traits
            i = i-1
    return worm_baby    

            