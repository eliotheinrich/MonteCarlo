import subprocess
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("job_name")
parser.add_argument("config")
parser.add_argument("--executable", nargs=1, default=["main"])
parser.add_argument("--mem", nargs=1, default=["5gb"])
parser.add_argument("--time", nargs=1, default=["24:00:00"])
parser.add_argument("--ncores", nargs=1, default=[64], type=int)
parser.add_argument("--partition", nargs=1, default=["default"])
args = parser.parse_args()

job_name = args.job_name
config_path = args.config
memory = args.mem[0]
time = args.time[0]
ncores = args.ncores[0]
executable = args.executable[0]
partition = args.partition[0]

if int(ncores) > 64:
	print("Cannot run a job with more than 64 CPU cores.")
	raise ValueError

if partition == "default":
	if int(ncores) > 48:
		partition = "full_nodes64"
	else:
		partition = "full_nodes48"



script = [f"#!/usr/bin/tcsh",
		  f"#SBATCH --partition={partition}",
		  f"#SBATCH --job-name={job_name}    # Job name",
		  f"#SBATCH  --ntasks 1 --cpus-per-task {ncores}",
		  f"#SBATCH --mem={memory}                     # Job memory request",
		  f"#SBATCH --time={time}               # Time limit hrs:min:sec",
		  f"cd ~/MonteCarlo/stable",
		  f"./{executable} /data/heinriea/MonteCarlo/{config_path} {ncores}"]

script = '\n'.join(script)

with open('job_tmp.sl', 'w') as f:
	f.write(script)

subprocess.run(["sbatch", "job_tmp.sl"])
subprocess.run(["rm", "job_tmp.sl"])
