#!/bin/bash
#SBATCH --job-name=quenching_cpp_actor
#SBATCH --account=hpc_p_spiteri
#SBATCH --nodes=2
#SBATCH --time=15-00:00:00
#SBATCH --mem=150G
#SBATCH --cpus-per-task=32
#SBATCH --constraint=cascade
#SBATCH --output=quench_new_actor_dynamic_test-%j.out

module --force purge
module load StdEnv/2020 gcc/9.3.0 openmpi/4.0.3
module load eigen/3.4.0 gsl/2.6 sundials/6.4.1

export CC=gcc
export CXX=g++
export LD_LIBRARY_PATH="/globalhome/tus210/HPC/lib64:$LD_LIBRARY_PATH"

echo "Building C++ CAF actor..."
make -f Makefile.cpp_dynamic clean all


if [ ! -f "./actor_cpp_bisection" ]; then
    echo "Error: actor_cpp_bisection not found after compilation"
    exit 1
fi
echo "Executable actor_cpp_bisection found"

nodes=( $(scontrol show hostname) )
server_node="${nodes[0]}"
client_nodes=( "${nodes[@]:1}" )

echo "Server node: $server_node"
echo "Client node(s): ${client_nodes[@]}"

##############################################################################

echo "Starting C++ CAF actor server on $server_node"
srun --nodes=1 --ntasks=1 --nodelist="${server_node}" ./actor_cpp_bisection -s -p 36494 --enable-early-termination --enable-bracket --caf.scheduler.max-threads=32 &

sleep 5  # Give server time to start

for client in "${client_nodes[@]}"; do
  echo "Starting client on $client"
  srun --nodes=1 --ntasks=1 --nodelist="$client" ./actor_cpp_bisection -p 36494 -H "$server_node"  --enable-early-termination --enable-bracket --caf.scheduler.max-threads=32 &
done

wait
echo "All processes have finished."
seff "$SLURM_JOB_ID"
