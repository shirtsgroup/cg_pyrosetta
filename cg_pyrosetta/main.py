import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(
        prog = 'CG_pyrosetta',
        description = 'Runs a simulated annealing MC simulation',
    )
    parser.add_argument(
        '--filename',
        action = 'store',
        default = 'input.params',
        metavar = 'f',
        desciption = 'the specific parameter file to use when running CG_pyrosetta',
    )
    parser.add_argument(


    )
