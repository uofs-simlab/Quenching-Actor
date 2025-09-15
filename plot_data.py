import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
from matplotlib.colors import LogNorm
from matplotlib.ticker import LogLocator, LogFormatter

def parse_out_file(file_path):
    data = []
    # Pattern to match lines like: "Received result for gg=3, theta=0, xs=140.08, uq=1.000"
    # or "Received result for gg=3, theta=0, xs=315.06, uq=-0.786461"
    old_pattern = re.compile(
        r"Received result for gg=(\d+), theta=([-?\d.]+), xs=([\d.]+), uq=([-?\d.eE]+)"
    )
    
    # Pattern to match lines with depth and neighbors (to avoid duplicates)
    # "Job completed: gg=0, theta=25.6129, xs=51.0572, result=-0.074517, runtime=541.84s, depth=3, neighbors=9"
    depth_pattern = re.compile(
        r"Job completed: gg=([^,]+), theta=([^,]+), xs=([^,]+), result=([^,]+), usmin=([^,]+), usmax=([^,]+), runtime=([^,]+)s, depth=([^,]+), neighbors=([^,]+)"
    )
    
    with open(file_path, 'r') as f:
        for line in f:
            # First try the depth pattern to avoid duplicates
            match = depth_pattern.search(line)
            if match:
                gg = int(float(match.group(1)))
                theta = float(match.group(2))
                xs = float(match.group(3))
                result_str = match.group(4)
                result_type = "result"
            else:
                # Try the old pattern (for older log files)
                match = old_pattern.search(line)
                if not match:
                    continue
                
                gg = int(match.group(1))
                theta = float(match.group(2))
                xs = float(match.group(3))
                result_str = match.group(4)
                result_type = "uq"
            
            try:
                result_val = float(result_str)
            except ValueError:
                result_val = np.nan
            
            # For visualization purposes:
            # - Treat result=1.0 or uq=1.0 (stable results) as NaN (don't plot them)
            # - Keep quenching results (result < 0 or uq < 0) as absolute values for log scale plotting
            if np.isclose(result_val, 1.0, atol=1e-6) or not np.isfinite(result_val):
                val = np.nan  # Don't plot stable results or invalid values
            elif result_val < 0:
                val = abs(result_val)  # Use absolute value for log scale (quenching results)
            else:
                val = np.nan  # Any other positive values that aren't 1.0, treat as NaN
            
            data.append((gg, theta, xs, val))
    
    return pd.DataFrame(data, columns=["gg", "theta", "xs", "result"])

def plot_heatmaps(df, output_dir="heatmaps", show_points=True): 
    os.makedirs(output_dir, exist_ok=True)
    unique_ggs = sorted(df["gg"].unique())
    gamma_map = {0: 0.025, 1: 0.02, 2: 0.01, 3: 0.001}

    if len(unique_ggs) == 0:
        print("No valid gg values found.")
        return

    uq_nonzero = df["result"][df["result"] > 0]
    if uq_nonzero.empty:
        print("All values are zero or NaN — cannot plot log scale.")
        return

    # Set fixed colorbar range from 0 to 10^3
    vmin = 1e-2  # Small positive value to avoid log(0)
    vmax = 1000   # 10^3 = 1000
    levels = np.logspace(np.log10(vmin), np.log10(vmax), 100)

    nrows, ncols = 2, 2
    fig, axes = plt.subplots(nrows, ncols, figsize=(10, 8), sharex=True, sharey=True, constrained_layout=True)
    axes = axes.flatten()

    contour_last = None

    for i, gg in enumerate(unique_ggs):
        ax = axes[i]
        full_subset = df[df["gg"] == gg]
        valid_subset = full_subset.dropna()
        nan_subset = full_subset[full_subset["result"].isna()]

        if valid_subset.empty:
            ax.set_visible(False)
            continue

        # Draw heatmap under everything
        contour = ax.tricontourf(
            valid_subset["xs"], valid_subset["theta"], valid_subset["result"],
            levels=levels, cmap="viridis", norm=LogNorm(vmin=vmin, vmax=vmax)
        )
        contour_last = contour

        # Plot points only if show_points is True
        if show_points:
            # Plot NaNs as very small gray crosses
            if not nan_subset.empty:
                ax.scatter(nan_subset["xs"], nan_subset["theta"], color='gray', marker='x', s=0.6)

            # Plot valid job points as very small black dots on top
            ax.scatter(valid_subset["xs"], valid_subset["theta"], color='black', marker='o', s=0.6)

        gamma_val = gamma_map.get(gg, f"gg={gg}")
        ax.set_title(f"$\\gamma = {gamma_val}$")


        if i // ncols == nrows - 1:
            ax.set_xlabel("$x_s$")
        if i % ncols == 0:
            ax.set_ylabel(r"$\theta$")

    for j in range(len(unique_ggs), len(axes)):
        axes[j].set_visible(False)

    if contour_last is not None:
        cbar = fig.colorbar(contour_last, ax=axes, orientation="vertical", shrink=0.8, aspect=30, pad=0.02)
        cbar.set_label("$|result|$")
        cbar.ax.tick_params(labelsize=10)
        cbar.ax.yaxis.set_major_locator(LogLocator(base=10.0, numticks=6))
        cbar.ax.yaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1, numticks=100))
        cbar.ax.yaxis.set_major_formatter(LogFormatter(labelOnlyBase=False))

    plt.savefig(f"{output_dir}/tricontourf_refinement2.png", dpi=300, bbox_inches='tight')
    plt.show()

def plot_single_gg(df, gg_value, output_dir="heatmaps", show_points=True, figsize=(10, 8)):
    """
    Plot a single gg value at a larger scale for detailed examination.
    
    Parameters:
    - df: DataFrame with the data
    - gg_value: The specific gg value to plot
    - output_dir: Directory to save the plot
    - show_points: Whether to show individual data points
    - figsize: Figure size tuple
    """
    os.makedirs(output_dir, exist_ok=True)
    gamma_map = {0: 0.025, 1: 0.02, 2: 0.01, 3: 0.001}
    
    # Filter data for the specific gg value
    subset = df[df["gg"] == gg_value]
    if subset.empty:
        print(f"No data found for gg={gg_value}")
        return
    
    valid_subset = subset.dropna()
    nan_subset = subset[subset["result"].isna()]
    
    if valid_subset.empty:
        print(f"No valid (quenching) data found for gg={gg_value}")
        # Still plot the NaN points if they exist
        if not nan_subset.empty:
            fig, ax = plt.subplots(figsize=figsize)
            ax.scatter(nan_subset["xs"], nan_subset["theta"], color='gray', marker='x', s=20, label='Not quenched')
            gamma_val = gamma_map.get(gg_value, f"gg={gg_value}")
            ax.set_title(f"$\\gamma = {gamma_val}$ (No quenching data)")
            ax.set_xlabel("$x_s$")
            ax.set_ylabel(r"$\theta$")
            ax.legend()
            plt.savefig(f"{output_dir}/single_gg_{gg_value}_no_quenching.png", dpi=300, bbox_inches='tight')
            plt.show()
        return
    
    # Set colorbar range
    vmin = max(1e-6, valid_subset["result"].min())
    vmax = min(1e3, valid_subset["result"].max())
    levels = np.logspace(np.log10(vmin), np.log10(vmax), 100)
    
    fig, ax = plt.subplots(figsize=figsize)
    
    # Create the contour plot
    contour = ax.tricontourf(
        valid_subset["xs"], valid_subset["theta"], valid_subset["result"],
        levels=levels, cmap="viridis", norm=LogNorm(vmin=vmin, vmax=vmax)
    )
    
    # Plot points if requested
    if show_points:
        # Plot NaN points as gray crosses
        if not nan_subset.empty:
            ax.scatter(nan_subset["xs"], nan_subset["theta"], color='gray', marker='x', s=15, alpha=0.7)
        
        # Plot valid points as black dots
        ax.scatter(valid_subset["xs"], valid_subset["theta"], color='black', marker='o', s=15, alpha=0.7)
    
    # Set title and labels
    gamma_val = gamma_map.get(gg_value, f"gg={gg_value}")
    ax.set_title(f"$\\gamma = {gamma_val}$ (gg={gg_value}) - Detailed View", fontsize=14)
    ax.set_xlabel("$x_s$", fontsize=12)
    ax.set_ylabel(r"$\theta$", fontsize=12)
    
    # Add colorbar
    cbar = fig.colorbar(contour, ax=ax, orientation="vertical", shrink=0.8, aspect=30)
    cbar.set_label("$|Uq|$", fontsize=12)
    cbar.ax.tick_params(labelsize=10)
    cbar.ax.yaxis.set_major_locator(LogLocator(base=10.0, numticks=6))
    cbar.ax.yaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1, numticks=100))
    cbar.ax.yaxis.set_major_formatter(LogFormatter(labelOnlyBase=False))
    
    # Print some statistics for this gg value
    print(f"\nStatistics for gg={gg_value}:")
    print(f"  Total points: {len(subset)}")
    print(f"  Quenching points: {len(valid_subset)}")
    print(f"  Stable points: {len(nan_subset)}")
    if len(valid_subset) > 0:
        print(f"  Quenching range: [{valid_subset['result'].min():.6f}, {valid_subset['result'].max():.6f}]")
        print(f"  xs range: [{valid_subset['xs'].min():.2f}, {valid_subset['xs'].max():.2f}]")
        print(f"  theta range: [{valid_subset['theta'].min():.2f}, {valid_subset['theta'].max():.2f}]")
    
    plt.tight_layout()
    plt.savefig(f"{output_dir}/single_gg_{gg_value}.png", dpi=300, bbox_inches='tight')
    plt.show()

if __name__ == "__main__":
    # file_path = "julia_quench_refine-3931789.out"
    file_path = "/globalhome/tus210/HPC/quenchin_actor/cpp_implementation/quench_cpp_static-4625162.out"  # Replace with your actual file path
    df = parse_out_file(file_path)
    
    # Print statistics
    print(f"Total data points: {len(df)}")
    num_nan = df["result"].isna().sum()
    num_quenching = len(df) - num_nan
    print(f"Number of stable results (NaN): {num_nan}")
    print(f"Number of quenching results: {num_quenching}")
    
    if num_quenching > 0:
        quenching_data = df.dropna()
        print(f"Quenching result range: [{quenching_data['result'].min():.6f}, {quenching_data['result'].max():.6f}]")
        
        # Show breakdown by gg value
        print("\nBreakdown by gg value:")
        for gg in sorted(df['gg'].unique()):
            gg_data = df[df['gg'] == gg]
            gg_quenching = gg_data.dropna()
            print(f"  gg={gg}: {len(gg_data)} total, {len(gg_quenching)} quenching")
    
    # Create plots
    plot_heatmaps(df)
    
    # Plot a single gg value at larger scale (example: gg=0)
    print("\n" + "="*50)
    print("Creating detailed plot for gg=0:")
    plot_single_gg(df, gg_value=0)
    
    # You can also create plots without showing the data points:
    # plot_heatmaps(df, show_points=False)
    # plot_single_gg(df, gg_value=0, show_points=False)
