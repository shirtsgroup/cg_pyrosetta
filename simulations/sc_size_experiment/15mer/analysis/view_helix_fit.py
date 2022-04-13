import dill
import matplotlib.pyplot as plt
import argparse


# argparse program that givne a directory iterates through interactive matplotlib sessions

def parse_args():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--files","-F",
        required = True,
        help = "list of files to view matplotlib pickled objects",
        nargs = '+'
    )
    
    return parser.parse_args()




def main():

    args = parse_args()

    for graph_file in args.files:
        print(graph_file)
        plot_info = dill.load(open(graph_file, 'rb'))

        fig = plt.figure()
        ax = fig.gca(projection='3d')
        # bb_helix
        ax.scatter3D(plot_info[0], plot_info[1], plot_info[2], s=125)
        ax.plot3D(plot_info[0], plot_info[1], plot_info[2])
        # fitted_helix
        ax.plot3D(plot_info[3], plot_info[4], plot_info[5])
        # center points + normal
        ax.scatter3D(plot_info[6], plot_info[7], plot_info[8], c="black")
        ax.plot3D(plot_info[9], plot_info[10], plot_info[11], c="black")
        # projected_points
        ax.scatter3D(plot_info[12], plot_info[13], plot_info[14], c="black")
        ax.plot3D(plot_info[12], plot_info[13], plot_info[14], "black")

        plt.show()





if __name__ == "__main__":
    main()
