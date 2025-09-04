#!/bin/bash
#SBATCH --job-name=quenching_exhaustive_early
#SBATCH --account=hpc_p_spiteri
#SBATCH --nodes=3
#SBATCH --time=5-00:00:00
#SBATCH --mem=50G
#SBATCH --cpus-per-task=32
#SBATCH --constraint=cascade
#SBATCH --output=julia_quench_exhustive_early-%j.out

module --force purge
module load StdEnv/2020 julia

export JULIA_PROJECT="/globalhome/tus210/HPC/quenchin_actor"
export LD_LIBRARY_PATH="/globalhome/tus210/HPC/lib64:$LD_LIBRARY_PATH"

echo "Setting up Julia environment on head node..."
julia --project="$JULIA_PROJECT" -e 'using Pkg; Pkg.resolve(); Pkg.instantiate()' 2>/dev/null || {
    echo "Julia setup failed!"
    exit 1
}

echo "Julia environment setup complete. Using exhaustive approach with early termination."

make -f Makefile.ex clean exhustive_test

if [ ! -f "./exhustive_test" ]; then
    echo "Error: exhustive_test not found after compilation"
    exit 1
fi
echo "Executable exhustive_test found"


nodes=( $(scontrol show hostname) )
server_node="${nodes[0]}"
client_nodes=( "${nodes[@]:1}" )

echo "Server node: $server_node"
echo "Client node(s): ${client_nodes[@]}"

srun --nodes=1 --ntasks=1 --nodelist="${server_node}" ./exhustive_test -s -p 2444 --caf.scheduler.max-threads=32 &

sleep 5  # Give server time to start

for client in "${client_nodes[@]}"; do
  echo "Starting client on $client"
  srun --nodes=1 --ntasks=1 --nodelist="$client" ./exhustive_test -p 2444 -H "$server_node" --caf.scheduler.max-threads=32 &
done

wait
echo "All processes have finished."

seff "$SLURM_JOB_ID"
