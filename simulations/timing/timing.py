import time
import subprocess
import numpy as np
import os

# Get general timing for running signac jobs

if os.path.isdir("workspace"):
    subprocess.run(["rm", "-rf", "workspace"])

# Output file

out_file = open("timing.txt", "w")

# Timing to initialize simulations

t1 = time.time()

subprocess.run(["python", "init.py"])

t2 = time.time()
delta_t = t2 - t1

print("Initialization took", delta_t, "seconds")
out_file.write("Initialization : " + str(delta_t) + " seconds\n")

# Timing changing parameters

t1 = time.time()

subprocess.run(["python", "project.py", "run", "-o", "set_parameters"])

t2 = time.time()
delta_t = t2 - t1

print("Changing parameters took", delta_t, "seconds")
out_file.write("Parameter Change : " + str(delta_t) + " seconds\n")

# Timing multiple simulation
times = []

for _ in range(25): # Should make this a variable that can change
    t1 = time.time()

    subprocess.run(["python", "project.py", "run", "-n", "1", "-o", "run_mc_simulation", "--debug"])

    t2 = time.time()
    delta_t = t2 - t1
    times.append(delta_t)
    print("3000 step simulation took", delta_t, "seconds")
    out_file.write("Simulation : " + str(delta_t) + " seconds\n")

print("Average 3000 step simulation took", np.mean(times), "seconds")
print("with a standard devation of", np.std(times))
print("Simulation runs", 3000/np.mean(times), "steps/second")
out_file.write("Average Simulation : " + str(delta_t) + " seconds\n")
out_file.write("Rate : " + str(3000/np.mean(times)) + " steps/seconds\n")

out_file.close()



