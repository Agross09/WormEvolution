import worm_gen_bit as worm_gen
import modified_bitsey_1 as bitsey
import sys



def main():
    trait_file = open(sys.argv[1], "r")
    population = worm_gen.gen_worms(trait_file)
    worm_to_sim = population.pop()
    print_worm(worm_to_sim)
    bitsey.setup_and_sim(worm_to_sim, 100)    
    trait_file.close()
    return

def print_worm(to_print):
    print("-----------------Worm to sim--------------------")
    print("Km = " + str(to_print.Km))
    print("N = " + str(to_print.N))
    print("GJ_scale = " + str(to_print.GJ_scale))
    print("num_cells = " + str(to_print.num_cells))
    print("G_k = " + str(to_print.G_k))
    print("G_k = " + str(to_print.G_na))
    print("G_cl = " + str(to_print.G_cl))
    print("GJ_diff_m = " + str(to_print.GJ_diff_m))
    print("\n")

main()
worm_evolution.py
Displaying main.py.