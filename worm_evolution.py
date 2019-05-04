import worm_gen 
import main as bitsey
import sys



def main():
    trait_file = open(sys.argv[1], "r")
    worm_gen.gen_worms(trait_file)
    
    trait_file.close()
    return

main()