import re
import matplotlib.pyplot as plt
import numpy as np

def extract_runtimes(filename):
    runtimes = []
    runtime_pattern = r'runtime=(\d+\.?\d*)s'
    
    try:
        with open(filename, 'r') as file:
            for line in file:
                match = re.search(runtime_pattern, line)
                if match:
                    runtime_value = float(match.group(1))
                    runtimes.append(runtime_value)
    except FileNotFoundError:
        print(f"Error: File {filename} not found.")
        return []
    
    return runtimes

def create_distribution_plot(dynamic_runtimes, static_runtimes):
    """Create a clean distribution comparison plot with average lines and counts."""
    
    plt.style.use('default')
    fig, ax = plt.subplots(figsize=(10, 6))
    
    max_runtime = max(max(dynamic_runtimes), max(static_runtimes))
    min_runtime = min(min(dynamic_runtimes), min(static_runtimes))
    bins = np.linspace(min_runtime, max_runtime, 60)
    
    plt.hist(dynamic_runtimes, bins=bins, alpha=0.7, label='Optimised dynamic regridding', 
             density=False, color='#87CEEB', edgecolor='none')
    plt.hist(static_runtimes, bins=bins, alpha=0.7, label='Static regridding', 
             density=False, color='#F08080', edgecolor='none')
    
    # Add average lines
    dynamic_average = np.mean(dynamic_runtimes)
    static_average = np.mean(static_runtimes)
    
    plt.axvline(dynamic_average, color='#4682B4', linestyle='--', linewidth=2, alpha=0.8,
                label=f'Dynamic mean: {dynamic_average:.1f}s')
    plt.axvline(static_average, color='#CD5C5C', linestyle='--', linewidth=2, alpha=0.8,
                label=f'baseline mean: {static_average:.1f}s')
    
    plt.yscale('log')

    plt.xlabel('Wall-clock time (seconds)', fontsize=12)
    plt.ylabel('Number of samples', fontsize=12)
    
    # Legend with average information
    plt.legend(fontsize=10, loc='upper right', frameon=True)

    plt.grid(False)
    
    # Set y-axis limits to be reasonable for log scale with counts
    plt.ylim(1, max(plt.ylim()[1], 10000))
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    plt.tight_layout()
    return fig

def main():
    # File paths
    dynamic_file = '/globalhome/tus210/HPC/quenchin_actor/cpp_implementation/quench_new_actor_dynamic_box_early_update-4662663.out'
    static_file = '/globalhome/tus210/HPC/quenchin_actor/cpp_implementation/quench_cpp_static_box-4645254.out'

    dynamic_runtimes = extract_runtimes(dynamic_file)
    static_runtimes = extract_runtimes(static_file)
    
    if not dynamic_runtimes or not static_runtimes:
        print("Could not extract runtime data from one or both files.")
        return
    
    print(f"Dynamic: {len(dynamic_runtimes)} samples")
    print(f"Static: {len(static_runtimes)} samples")
    
    fig = create_distribution_plot(dynamic_runtimes, static_runtimes)
    
    output_filename = 'runtime_distributions_logscale.pdf'
    plt.savefig(output_filename, format='pdf', bbox_inches='tight', facecolor='white', edgecolor='none')
    print(f"Log-scale plot saved as: {output_filename}")
    
    # Print basic statistics including averages
    print(f"\nDynamic Implementation:")
    print(f"  Mean: {np.mean(dynamic_runtimes):.2f}s, Median: {np.median(dynamic_runtimes):.2f}s")
    print(f"  Std: {np.std(dynamic_runtimes):.2f}s")
    print(f"Static Implementation:")
    print(f"  Mean: {np.mean(static_runtimes):.2f}s, Median: {np.median(static_runtimes):.2f}s")
    print(f"  Std: {np.std(static_runtimes):.2f}s")
    
    plt.show()

if __name__ == "__main__":
    main()