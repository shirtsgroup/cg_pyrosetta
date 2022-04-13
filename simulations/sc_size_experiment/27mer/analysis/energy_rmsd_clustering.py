import argparse
import signac
import cg_pyrosetta
import numpy as np
import pandas as pd
import copy
import os
import sys
import shutil
import flow
import analyze_foldamers
import matplotlib.pyplot as plt
from matplotlib import cm
import mdtraj as md
from flow import FlowProject

plt.rcParams["font.family"] = "serif"
plt.rcParams.update({"font.size": 15})


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--signac_dir",
        help="Directory containing signac project",
        required=False,
        type=str,
        default=os.path.abspath(""),
    )
    parser.add_argument(
        "--output_dir",
        help="Output directory for generating figures",
        required=False,
        type=str,
        default=os.path.abspath("output"),
    )
    parser.add_argument(
        "--rerun",
        action="store_true",
        help="Flag to rerun analysis on parameter sets if files directory already exists",
    )
    parser.add_argument(
        "--single",
        help="Run analysis on as single set of parameters",
        type=str,
        default=None,
    )
    parser.add_argument(
        "--cluster_parameters",
        help="specify key-value pairs of which parameters to feed the "
        + "analyze_foldamer.get_cluster_medoid_positions_DBSCAN()"
        + "function\n"
        + "Possible keys are: eps, min_samples, frame_start, frame_stride,"
        + " frame_end, output_format, ride, frame_end, output_format, "
        + "cgmodel, plot_silhouette, filter, filter_ratio, "
        + "output_cluster_traj, core_points_only, np",
        nargs="+",
        default=None,
    )
        if os.path.isfile("cluster_parameters.yml"):

    args = parser.parse_args()
    return args


def get_commandline_input():
    cmd_line = sys.argv
    cmd_line.insert(0, sys.executable)
    cmd_str = " ".join(cmd_line)
    return cmd_str


def get_default_parameters():
    default_params = {
        "eps": 0.15,
        "min_samples": 30,
        "frame_start": 700,
        "frame_stride": 1,
        "frame_end": -1,
        "output_format": "pdb",
        "plot_silhouette": True,
        "filter": True,
        "filter_ratio": 0.5,
        "output_cluster_traj": True,
        "core_points_only": False,
        "parallel" : -1,
    }
    return default_params

def compute_2_population_z_score(X1, X2):
    mu_1 = np.mean(X1)
    mu_2 = np.mean(X2)
    sigma_1 = np.std(X1)
    sigma_2 = np.std(X2)

    z_score = (mu_1 - mu_2) / np.sqrt(np.square(sigma_1)/len(X1) + np.square(sigma_2)/len(X2)) 
    return z_score


def main():
    args = parse_args()
    os.chdir(args.signac_dir)

    # Getting Signac project schema
    print("Fetching signac schema information...")
    project = signac.get_project()
    schema = project.detect_schema()
    sc_sizes, reps = schema.items()
    sc_sizes = list(sc_sizes[1][float])
    sc_sizes.sort()
    select_sc_sizes = [str(round(a, 3)) for a in sc_sizes]
    print("List of all side chain sizes:", *select_sc_sizes)

    # Get commandline scirpt
    cmd_str = get_commandline_input()

    # Get parameters for clustering (Update from default values)
    param_dict = get_default_parameters()
    if args.cluster_parameters is not None:
        if len(args.cluster_parameters) % 2 != 0:
            print(
                "One unmatched key-value pairs. Remvoing",
                args.cluster_parameters.pop(),
                "from input parameters.",
            )
            print("Remaining input key-value pairs are:", *args.cluster_parameters)

        for key, value in zip(
            args.cluster_parameters[::2], args.cluster_parameters[1::2]
        ):
            if key in param_dict.keys():
                param_dict[key] = type(param_dict[key])(value)
        print("Using non-default parameters for clutsering:")
        print(param_dict)

    # Single parameter option
    if args.single is not None:
        i_single = select_sc_sizes.index(args.single)
        sc_sizes = [sc_sizes[i_single]]

    # Setup output directory
    if not os.path.isdir(args.output_dir):
        os.mkdir(args.output_dir)

    # kT values for tempeature plot
    kts = np.array([10 * ( 0.9 ) ** i for i in range(50)])
    out_steps = np.array([10000 * i for i in range(50)])
    
    for sc_size in sc_sizes:
        print("Working on SC Size =", str(round(sc_size, 3)), "...")
        # Setup parameters for each set of unique simulations
        sub_dir = os.path.join(args.output_dir, "sc_size_" + str(round(sc_size, 3)))
        if os.path.isdir(sub_dir):
            if args.rerun:
                shutil.rmtree(sub_dir)
            else:
                print(
                    "Skipping",
                    sub_dir,
                    "if you want to rerun this analysis add the --rerun flag",
                )
                continue
        os.mkdir(sub_dir)

        # Command line that generated specific output
        with open(os.path.join(sub_dir, "command_line_call.txt"), "w") as fw:
            fw.write(cmd_str + "\n")

        # Get trajectory and energies of each simulation
        traj_file_list = []
        all_energies = []
        energy_traj = []
        for job in project.find_jobs({"sc_size": sc_size}):
            structure_file = job.fn("trajectory.pdb")
            traj_file_list.append(structure_file)
            energy_file = job.fn("energies.txt")
            energies = pd.read_csv(energy_file, header=None)
            energy_traj.append(energies.values[:, 1])
            all_energies.extend(energies.values[:, 1])
        all_energies = np.array(all_energies)
        
        # Energy trajectory plot
        fig, ax1 = plt.subplots(1,1,figsize = [20,10])
        ax2 = ax1.twinx()
        ax2.plot(out_steps, kts, 'r')
        ax2.set_xlim([0, 500000])
        ax2.set_ylabel("Simulated Temperature", color = 'r')
        for traj in energy_traj:
            ax1.plot(energies.values[:,0], traj, alpha = 0.4, lw=2)
            ax1.set_xlabel("Steps")
            ax1.set_ylabel("Energy (A.E.U)")
            ax1.set_title("SC size = " + str(round(job.sp['sc_size'],4)))

        fig.savefig(os.path.join(sub_dir, "energy_trajectory.jpg"), bbox_inches="tight")
        fig.savefig(os.path.join(sub_dir, "energy_trajectory.pdf"), bbox_inches="tight")

        plt.close("all")

        # Run clustering
        (
            medoid_positions,
            cluster_sizes,
            cluster_rmsd,
            n_noise,
            silhouette_avg,
            labels,
            original_indices,
        ) = analyze_foldamers.cluster.get_cluster_medoid_positions_DBSCAN(
            traj_file_list,
            None,
            output_dir=os.path.join(sub_dir, "cluster_output"),
            **param_dict
        )

        if len(cluster_sizes) > 0:
            # Output cluster energy distributions
            clusters = list(np.unique(labels))
            all_cluster_energies = []
            if -1 in clusters:
                clusters.remove(-1)
            for i in clusters:
                plt.figure(figsize=[5, 5], dpi=100)
                cluster_indices = np.where(labels == i)[0]
                cluster_energies = all_energies[original_indices[cluster_indices]]
                all_cluster_energies.append(cluster_energies)
                mean_energy = np.mean(cluster_energies)
                std_energy = np.std(cluster_energies)
                color = cm.nipy_spectral(float(i) / len(clusters))
                plt.hist(cluster_energies, bins=50, color="red")
                plt.xlabel("Mean Energy (A.E.U)")
                plt.ylabel("Energy Std. Err. (A.E.U)")
                plt.title("Cluster " + str(i) + " Energy Distribution")
                plt.savefig(
                    os.path.join(
                        sub_dir, "cluster_energy_dist_cluster_" + str(i) + ".pdf"
                    ), bbox_inches="tight"
                )
                plt.savefig(
                    os.path.join(
                        sub_dir, "cluster_energy_dist_cluster_" + str(i) + ".jpg"
                    ), bbox_inches="tight"
                )
                plt.close("all")

            fig = plt.figure(figsize=[10, 10], dpi=500)
            ax = fig.add_axes([0, 0, 1, 1])
            energy_dists = []
            for i in clusters:
                cluster_indices = np.where(labels == i)[0]
                cluster_energies = all_energies[original_indices[cluster_indices]]
                energy_dists.append(cluster_energies)

            parts = ax.violinplot(
                energy_dists,
                positions=cluster_rmsd,
                showextrema=False,
                showmeans=False,
                widths=0.02,
            )
            for i, pc in enumerate(parts["bodies"]):
                color = cm.nipy_spectral(float(i) / len(clusters))
                pc.set_facecolor(color)
                pc.set_edgecolor("black")
                pc.set_alpha(0.5)
            # quartile1, medians, quartile3 = np.percentile(energy_dists, [25, 50, 75])
            # plt.scatter(medians, marker = 'o')
            plt.xlabel("Cluster RMSF (A.L.U)")
            plt.ylabel("Cluster Mean Energy (A.E.U)")
            plt.legend(["Cluster " + str(a) for a in clusters], loc="best")
            plt.savefig(os.path.join(sub_dir, "rmsd_energy_violin_plot.pdf"), bbox_inches="tight")
            plt.savefig(os.path.join(sub_dir, "rmsd_energy_violin_plot.jpg"), bbox_inches="tight")
            
            plt.close("all")

            # RMSD scatter plot vs energies

            # Getting a trajectory with all frames
            sim_traj = []
            all_traj = None
            for job in project.find_jobs({'sc_size':sc_size}):
                job_traj = md.load(job.fn("trajectory.pdb"))
                sim_traj.append(job_traj)
                if all_traj is None:
                    all_traj = job_traj
                else:
                    all_traj = all_traj.join(job_traj)

            plt.figure(figsize = [10, 10], dpi=500)
            energy_stdevs = []
            for i in clusters:
                medoid_mdtraj = md.load(os.path.join(sub_dir, "cluster_output/medoid_" + str(i) + ".pdb"))
                cluster_indices = np.where(labels == i)[0]
                rmsds_medoid_i = md.rmsd(all_traj[original_indices[cluster_indices]], medoid_mdtraj)
                cluster_energies = all_energies[original_indices[cluster_indices]]
                color = cm.nipy_spectral(float(i) / len(clusters))
                plt.scatter(rmsds_medoid_i, cluster_energies, s=1, alpha=0.7, c=color)
                plt.xlabel("RMSD to Cluster Medoid (A.L.U)")
                plt.ylabel("Cluster Energies (A.E.U)")
            plt.legend(["Cluster "+ str(a) for a in clusters], loc = "best")
            plt.savefig(os.path.join(sub_dir, "rmsd_energy_scatter_plot.pdf"), bbox_inches="tight")
            plt.savefig(os.path.join(sub_dir, "rmsd_energy_scatter_plot.jpg"), bbox_inches="tight")

            plt.close("all")

            # Distance matrix of each medoid structure to one another

            medoid_rmsd_matrix = 1

            cluster_output_files = os.listdir(os.path.join(sub_dir, "cluster_output"))
            cluster_output_files.sort()
            medoid_traj = None
            for medoid_file in cluster_output_files:
                if "medoid" in medoid_file:
                    if medoid_traj is None:
                        medoid_traj = md.load(os.path.join(sub_dir, "cluster_output", medoid_file))
                    else:
                        medoid_traj = medoid_traj.join(md.load(os.path.join(sub_dir, "cluster_output", medoid_file)))

            medoid_traj.superpose(medoid_traj)
            medoid_rmsd_matrix = np.zeros((medoid_traj.n_frames, medoid_traj.n_frames))

            for i in range(medoid_traj.n_frames):
                medoid_rmsd_matrix[i, :] = md.rmsd(medoid_traj, medoid_traj[i])

            plt.figure(figsize = [10,10])
            plt.matshow(medoid_rmsd_matrix, cmap = "viridis")
            ax = plt.gca()

            labels = ["Medoid " + str(i) for i in clusters]
            labels.insert(0, "")
            ax.set_xticklabels(labels)
            ax.set_yticklabels(labels)

            # for i in range(len(clusters)):
            #    for j in range(len(clusters)):
            #        plt.text(j, i, round(medoid_rmsd_matrix[i, j], 3), horizontalalignment='center', verticalalignment='center')


            plt.savefig(os.path.join(sub_dir, "intermedoid_rmsd_matrix.pdf"), bbox_inches="tight")
            plt.savefig(os.path.join(sub_dir, "intermedoid_rmsd_matrix.pdf"), bbox_inches="tight")
            plt.close("all")

            print(["Medoid " + str(i) for i in clusters])

            # Energy statistical test 
            # Z test for comparing two population means

            z_score_matrix = np.zeros((len(clusters), len(clusters)))
            for i in clusters:
                for j in clusters:
                    if i != j:
                        z_score_matrix[i,j] = compute_2_population_z_score(all_cluster_energies[i], all_cluster_energies[j])

            np.savetxt(os.path.join(sub_dir,"energy_z_scores.txt"), z_score_matrix, delimiter=',')
            print(z_score_matrix)

        else:
            print("No cluster Identified")


if __name__ == "__main__":
    main()
