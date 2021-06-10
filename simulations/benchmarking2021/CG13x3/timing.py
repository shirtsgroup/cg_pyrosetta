import time
import subprocess
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

# Timing single simulation

t1 = time.time()

subprocess.run(["python", "project.py", "run", "-n", "1", "-o", "run_mc_simulation"])

t2 = time.time()
delta_t = t2 - t1

print("A 3000 step simulation took", delta_t, "seconds")
out_file.write("Simulation : " + str(delta_t) + " seconds\n")

print("Simulation runs", 3000/delta_t, "steps/second")
out_file.write("Rate : " + str(3000/delta_t) + " steps/seconds\n")

out_file.close()



