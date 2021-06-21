import signac
import cg_pyrosetta
import os
import argparse
import warnings
import shutil as sh

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--statepoint",
        help="Define which statepoints to restart. State points can be defines the folowing ways:\n"
        + "1. single values"
        + "2. ranges denoted by - (100-200)"
        + "3. ranges denoted by <, >, (<500)",
        nargs="+",
        required = True,
        default=None,
    )

    args = parser.parse_args()
    return args


def process_arguments(statepoint_list, statepoints):
    """
    """

    # Get indexes of statepoint in input string list
    if len(statepoint_list) == 0:
        warnings.warn("Warning: No provided state points were "
        + "defined in the current signac projec\n",
        + "Available statepoints are: " + ", ".join(statepoints))

    if len(statepoint_list) % 2 == 1:
        warnings.warn("Warning: Only single values/logicals for"
        + "statepoint definitions\n",
        + "Be sure you enter a single entry for each statepoint")

    # Make selection dictionary
    sp_dict = {}
    for i in range(int(len(statepoint_list)/2)):
        if statepoint_list[2*i] in statepoints:
            if ">" not in statepoint_list[2*i+1] and \
            "<" not in statepoint_list[2*i+1]:
                sp_dict[statepoint_list[2*i]] = float(statepoint_list[2*i+1])
            else:
                if ">" in statepoint_list[2*i+1]:
                    sp_dict[statepoint_list[2*i]+".$gt"] = float(statepoint_list[2*i+1][1:])
                if "<" in statepoint_list[2*i+1]:
                    sp_dict[statepoint_list[2*i]+".$lt"] = float(statepoint_list[2*i+1][1:])
        else:
            warnings.warn("Warning: " + statepoint_list[2*i] + " is not "
                + "a valid state point defined in the current signac "
                + "projec\n"
                + "Available statepoints are: " + ", ".join(statepoints))

    return sp_dict

def main():
    # Parse arguments
    args = parse_args()
    
    # Get statepoints from signac
    project = signac.get_project()
    project_statepoints = [*project.detect_schema()]
    selection = process_arguments(args.statepoint, project_statepoints)

    print(selection)
    for job in project.find_jobs(selection):
        print(job.sp)
        job.remove()
        job.init()
        sh.copytree("parameters", os.path.join(job.workspace(), "parameters"))



if __name__ == "__main__":
    main()