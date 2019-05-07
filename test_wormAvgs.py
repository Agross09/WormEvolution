
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

