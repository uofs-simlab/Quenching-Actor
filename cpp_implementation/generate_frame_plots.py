#!/usr/bin/env python3
"""
Simple visualization script for FHN snapshots
Creates individual PNG files for each time frame, showing only u variable
Organized in folders by perturbation value
"""

import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import re

def load_snapshot_data(filename):
    """Load snapshot data from file"""
    data = []
    with open(filename, 'r') as f:
        time_line = f.readline().strip()
        time = float(time_line.split(':')[1])
        
        f.readline()  # Skip column header
        
        for line in f:
            if line.strip():
                parts = line.strip().split()
                xi, u, v = float(parts[0]), float(parts[1]), float(parts[2])
                data.append([xi, u, v])
    
    data = np.array(data)
    return time, data[:, 0], data[:, 1], data[:, 2]  # time, xi, u, v

def process_perturbation(perturbation_value):
    """Process all snapshots for a given perturbation value"""
    
    # Create output directory for this perturbation (new parameter set)
    if perturbation_value >= 0:
        folder_name = f"plots_xs5_theta10_A_{perturbation_value:.2f}".replace('.', 'p')
        pattern = f"fhn_snapshots_A{perturbation_value:.4f}*".replace('.', 'p')
    else:
        abs_val = abs(perturbation_value)
        folder_name = f"plots_xs5_theta10_A_n{abs_val:.2f}".replace('.', 'p')
        pattern = f"fhn_snapshots_An{abs_val:.3f}*".replace('.', 'p')
    
    os.makedirs(folder_name, exist_ok=True)
    
    print(f"Processing perturbation A = {perturbation_value}")
    print(f"Looking for pattern: {pattern}")
    
    # Find all snapshot files for this perturbation
    files = glob.glob(pattern)
    files.sort(key=lambda x: int(re.findall(r'frame_(\d+)', x)[0]))
    
    if not files:
        print(f"No files found for perturbation {perturbation_value}")
        return
    
    print(f"Found {len(files)} snapshot files")
    print(f"Output folder: {folder_name}")
    
    # Process each snapshot
    for i, filename in enumerate(files):
        if i % 20 == 0:
            print(f"Processing file {i+1}/{len(files)}: {filename}")
        
        time, xi, u, v = load_snapshot_data(filename)
        
        # Create figure for u component only
        plt.figure(figsize=(12, 6))
        plt.plot(xi, u, 'b-', linewidth=1.5, label='u')
        plt.title(f'FHN u component (xs=5, θ=10, A = {perturbation_value}, t = {time:.3f})', fontsize=14)
        plt.xlabel('ξ', fontsize=12)
        plt.ylabel('u', fontsize=12)
        plt.grid(True, alpha=0.3)
        
        # Set consistent y-axis limits for comparison across frames
        plt.ylim(-2, 2)  # Adjust as needed based on your data range
        
        # Save the figure
        frame_number = int(re.findall(r'frame_(\d+)', filename)[0])
        output_filename = f"{folder_name}/frame_{frame_number:03d}_t_{time:.3f}.png"
        plt.savefig(output_filename, dpi=150, bbox_inches='tight')
        plt.close()
    
    print(f"Saved {len(files)} individual frames to {folder_name}/")

def main():
    """Main function to generate all plots"""
    
    print("=== FHN Individual Frame Generator ===")
    
    # Perturbations to analyze
    perturbations = [0.0, -0.58, -0.5]
    
    # Process each perturbation
    for A in perturbations:
        process_perturbation(A)
        print()  # Empty line between perturbations
    
    print("=== Frame generation complete ===")
    
    # Summary
    for A in perturbations:
        if A >= 0:
            folder_name = f"plots_xs5_theta10_A_{A:.2f}".replace('.', 'p')
        else:
            abs_val = abs(A)
            folder_name = f"plots_xs5_theta10_A_n{abs_val:.2f}".replace('.', 'p')
        
        if os.path.exists(folder_name):
            num_files = len([f for f in os.listdir(folder_name) if f.endswith('.png')])
            print(f"  {folder_name}: {num_files} PNG files")

if __name__ == "__main__":
    main()