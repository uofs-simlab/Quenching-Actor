#!/usr/bin/env python3
"""
Simple Runtime Analysis Script
Shows only file comparison histograms for quenching vs non-quenching cases
"""

import re
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def parse_runtime_data(file_path):
    """
    Parse runtime data from output file
    Returns lists of runtimes for quenching and non-quenching cases
    """
    quenching_runtimes = []
    non_quenching_runtimes = []
    
    with open(file_path, 'r') as f:
        for line in f:
            # Look for job completion lines with runtime and uq values
            match = re.search(r'Job completed:.*uq=([+-]?\d+\.?\d*),.*runtime=(\d+\.?\d*)s', line)
            if match:
                uq_value = float(match.group(1))
                runtime = float(match.group(2))
                
                # uq=1.0 means no quenching, negative values mean quenching occurred
                if abs(uq_value - 1.0) < 1e-6:  # Non-quenching case
                    non_quenching_runtimes.append(runtime)
                else:  # Quenching case
                    quenching_runtimes.append(runtime)
    
    return quenching_runtimes, non_quenching_runtimes

def analyze_and_plot():
    """
    Simple analysis with only the requested histograms
    """
    # File paths
    dynamic_file = '/globalhome/tus210/HPC/quenchin_actor/cpp_implementation/quench_new_actor_dynamic_box_early_update-4662663.out'
    static_file = '/globalhome/tus210/HPC/quenchin_actor/cpp_implementation/quench_cpp_static_box-4645254.out'
    
    print("SIMPLE RUNTIME ANALYSIS")
    print("="*40)
    
    # Parse data from both files
    print(f"Processing files...")
    dynamic_quench, dynamic_no_quench = parse_runtime_data(dynamic_file)
    static_quench, static_no_quench = parse_runtime_data(static_file)
    
    # Print simple statistics
    print(f"\nDynamic file - Non-quenching: {len(dynamic_no_quench)} cases, Mean: {np.mean(dynamic_no_quench):.2f}s")
    print(f"Dynamic file - Quenching: {len(dynamic_quench)} cases, Mean: {np.mean(dynamic_quench):.2f}s")
    print(f"Static file - Non-quenching: {len(static_no_quench)} cases, Mean: {np.mean(static_no_quench):.2f}s")
    print(f"Static file - Quenching: {len(static_quench)} cases, Mean: {np.mean(static_quench):.2f}s")
    
    # Create only the two requested plots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle('Runtime Comparison: Dynamic vs Static Implementation', fontsize=14)
    
    alpha = 0.7
    
    # Non-quenching cases comparison (Dynamic vs Static)
    ax1.hist(dynamic_no_quench, bins=30, alpha=alpha, color='#3498db',
             label=f'Dynamic (n={len(dynamic_no_quench)})')
    ax1.hist(static_no_quench, bins=30, alpha=alpha, color='#e74c3c',
             label=f'Static (n={len(static_no_quench)})')
    ax1.axvline(np.mean(dynamic_no_quench), color='#3498db', linestyle='--', linewidth=2,
                label=f'Dynamic Mean: {np.mean(dynamic_no_quench):.1f}s')
    ax1.axvline(np.mean(static_no_quench), color='#e74c3c', linestyle='--', linewidth=2,
                label=f'Static Mean: {np.mean(static_no_quench):.1f}s')
    ax1.set_title('Non-Quenching Cases')
    ax1.set_xlabel('Runtime (seconds)')
    ax1.set_ylabel('Frequency')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Quenching cases comparison (Dynamic vs Static)
    ax2.hist(dynamic_quench, bins=30, alpha=alpha, color='#3498db',
             label=f'Dynamic (n={len(dynamic_quench)})')
    ax2.hist(static_quench, bins=30, alpha=alpha, color='#e74c3c',
             label=f'Static (n={len(static_quench)})')
    ax2.axvline(np.mean(dynamic_quench), color='#3498db', linestyle='--', linewidth=2,
                label=f'Dynamic Mean: {np.mean(dynamic_quench):.1f}s')
    ax2.axvline(np.mean(static_quench), color='#e74c3c', linestyle='--', linewidth=2,
                label=f'Static Mean: {np.mean(static_quench):.1f}s')
    ax2.set_title('Quenching Cases')
    ax2.set_xlabel('Runtime (seconds)')
    ax2.set_ylabel('Frequency')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # Save the plot
    output_file = 'simple_runtime_histograms.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\nHistograms saved as: {output_file}")
    
    plt.show()

if __name__ == "__main__":
    analyze_and_plot()