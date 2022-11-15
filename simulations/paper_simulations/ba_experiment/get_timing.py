import signac
import os
import numpy as np



def main():
    project = signac.get_project()
    run_times = []
    for job in project.find_jobs():
        if os.path.isfile(job.fn("minimum.pdb")):
            run_times.append(job.document.run_time)
    
    print("Runtime:", round(np.mean(run_times)/60,3), "m")
    print("Runtime STD:", round(np.std(run_times)/60,3), "m")

if __name__ == "__main__":
    main()