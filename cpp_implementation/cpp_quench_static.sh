#!/bin/bash
#SBATCH --job-name=quenching_cpp_static
#SBATCH --account=hpc_p_spiteri
#SBATCH --nodes=2
#SBATCH --time=10-00:00:00
#SBATCH --mem=5G
#SBATCH --cpus-per-task=32
#SBATCH --constraint=cascade
#SBATCH --output=quench_cpp_static-%j.out

module --force purge
module load StdEnv/2020 gcc/9.3.0 openmpi/4.0.3
module load eigen/3.4.0 gsl/2.6 sundials/6.4.1

export CC=gcc
export CXX=g++
export LD_LIBRARY_PATH="/globalhome/tus210/HPC/lib64:$LD_LIBRARY_PATH"

echo "Building C++ CAF static actor using Makefile..."
make -f Makefile.cpp_static

if [ $? -ne 0 ]; then
    echo "Compilation failed!"
    exit 1
fi

if [ ! -f "./actor_cpp_static" ]; then
    echo "Error: actor_cpp_static not found after compilation"
    exit 1
fi
echo "Executable actor_cpp_static found"

echo "Setting runtime library path..."
export LD_LIBRARY_PATH="/globalhome/tus210/HPC/lib64:$LD_LIBRARY_PATH"
echo "LD_LIBRARY_PATH: $LD_LIBRARY_PATH"


nodes=( $(scontrol show hostname) )
server_node="${nodes[0]}"
client_nodes=( "${nodes[@]:1}" )

echo "Server node: $server_node"
echo "Client node(s): ${client_nodes[@]}"

echo "Starting C++ CAF static actor server on $server_node"
srun --nodes=1 --ntasks=1 --nodelist="${server_node}" --export=ALL ./actor_cpp_static -s -p 32444 --caf.scheduler.max-threads=32 &

sleep 5  # Give server time to start

for client in "${client_nodes[@]}"; do
  echo "Starting client on $client"
  srun --nodes=1 --ntasks=1 --nodelist="$client" --export=ALL ./actor_cpp_static -p 32444 -H "$server_node" --caf.scheduler.max-threads=32 &
done

wait
echo "All processes have finished."

seff "$SLURM_JOB_ID"
