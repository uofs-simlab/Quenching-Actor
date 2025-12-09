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
    old_pattern = re.compile(
        r"Received result for gg=(\d+), theta=([-?\d.]+), xs=([\d.]+), uq=([-?\d.eE]+)"
    )
    
    depth_pattern = re.compile(
        r"Job completed: gg=([^,]+), theta=([^,]+), xs=([^,]+), uq=([^,]+), bracket=\[([^\]]+)\], runtime=([^,]+)s"
    )
    
    with open(file_path, 'r') as f:
        for line in f:
            match = depth_pattern.search(line)
            if match:
                gg = int(float(match.group(1)))
                theta = float(match.group(2))
                xs = float(match.group(3))
                result_str = match.group(4)
                result_type = "result"
            else:
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
                val = np.nan  
            elif result_val < 0:
                val = abs(result_val)  
            else:
                val = np.nan 
            
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
    fig, axes = plt.subplots(nrows, ncols, figsize=(10, 8), sharex=False, sharey=True, constrained_layout=True)
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

        contour = ax.tricontourf(
            valid_subset["xs"], valid_subset["theta"], valid_subset["result"],
            levels=levels, cmap="viridis", norm=LogNorm(vmin=vmin, vmax=vmax)
        )
        contour_last = contour

        if show_points:
            if not nan_subset.empty:
                ax.scatter(nan_subset["xs"], nan_subset["theta"], color='gray', marker='x', s=0.2)

            ax.scatter(valid_subset["xs"], valid_subset["theta"], color='black', marker='o', s=0.2)

        gamma_val = gamma_map.get(gg, f"gg={gg}")
        ax.set_title(f"$\\gamma = {gamma_val}$")
        
        # Set x-axis to start from exactly 0
        ax.set_xlim(left=-10.0)

        if i // ncols == nrows - 1:
            ax.set_xlabel(r"$x_s$")
        if i % ncols == 0:
            ax.set_ylabel(r"$\theta$")

    for j in range(len(unique_ggs), len(axes)):
        axes[j].set_visible(False)
    
    # Set x-axis limits for all visible plots to start from 0
    for i in range(len(unique_ggs)):
        axes[i].set_xlim(left=-10.0)

    if contour_last is not None:
        cbar = fig.colorbar(contour_last, ax=axes, orientation="vertical", shrink=0.8, aspect=30, pad=0.02)
        cbar.set_label(r"$|A^{\star}|$")
        cbar.ax.tick_params(labelsize=10)
        cbar.ax.yaxis.set_major_locator(LogLocator(base=10.0, numticks=6))
        cbar.ax.yaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1, numticks=100))
        cbar.ax.yaxis.set_major_formatter(LogFormatter(labelOnlyBase=False))

    plt.savefig(f"{output_dir}/tricontourf_refinement.pdf", dpi=300, bbox_inches='tight')
    plt.show()

# def plot_all_separate_figures(df, output_dir="heatmaps", show_points=True, figsize=(8, 6)):

#     os.makedirs(output_dir, exist_ok=True)
#     unique_ggs = sorted(df["gg"].unique())
    
#     if len(unique_ggs) == 0:
#         print("No valid gg values found.")
#         return
    
    # Plot each gg value separately
    # for gg_value in unique_ggs:
    #     print(f"\nCreating separate plot for gg={gg_value}:")
    #     plot_single_gg(df, gg_value=gg_value, output_dir=output_dir, 
    #                   show_points=show_points, figsize=figsize)

# def plot_single_gg(df, gg_value, output_dir="heatmaps", show_points=True, figsize=(8, 6)):
 
#     os.makedirs(output_dir, exist_ok=True)
#     gamma_map = {0: 0.025, 1: 0.02, 2: 0.01, 3: 0.001}
    
#     subset = df[df["gg"] == gg_value]
#     if subset.empty:
#         print(f"No data found for gg={gg_value}")
#         return
    
#     valid_subset = subset.dropna()
#     nan_subset = subset[subset["result"].isna()]
    
#     if valid_subset.empty:
#         print(f"No valid (quenching) data found for gg={gg_value}")
#         if not nan_subset.empty:
#             fig, ax = plt.subplots(figsize=figsize)
#             ax.scatter(nan_subset["xs"], nan_subset["theta"], color='gray', marker='x', s=0.5, label='Not quenched')
#             gamma_val = gamma_map.get(gg_value, f"gg={gg_value}")
#             ax.set_title(f"$\\gamma = {gamma_val}$ (No quenching data)")
#             ax.set_xlabel("$x_s$")
#             ax.set_ylabel(r"$\theta$")
#             ax.legend()
#             plt.savefig(f"{output_dir}/single_gg_{gg_value}_no_quenching.pdf", dpi=300, bbox_inches='tight')
#             plt.show()
#         return
    
#     # Set fixed colorbar range from 0 to 10^3 (consistent with combined plot)
#     vmin = 1e-2  # Small positive value to avoid log(0)
#     vmax = 1000   # 10^3 = 1000
#     levels = np.logspace(np.log10(vmin), np.log10(vmax), 100)
    
#     fig, ax = plt.subplots(figsize=figsize)
    
#     contour = ax.tricontourf(
#         valid_subset["xs"], valid_subset["theta"], valid_subset["result"],
#         levels=levels, cmap="viridis", norm=LogNorm(vmin=vmin, vmax=vmax)
#     )
    
#     if show_points:
#         if not nan_subset.empty:
#             ax.scatter(nan_subset["xs"], nan_subset["theta"], color='gray', marker='x', s=0.2, alpha=0.7)
        
#         ax.scatter(valid_subset["xs"], valid_subset["theta"], color='black', marker='o', s=0.2, alpha=0.7)
    
#     gamma_val = gamma_map.get(gg_value, f"gg={gg_value}")
#     ax.set_title(f"$\\gamma = {gamma_val}$ (gg={gg_value}) - Detailed View", fontsize=14)
#     ax.set_xlabel("$x_s$", fontsize=12)
#     ax.set_ylabel(r"$\theta$", fontsize=12)
    
#     ax.set_xlim(left=-10.0)
    
#     # Add colorbar
#     cbar = fig.colorbar(contour, ax=ax, orientation="vertical", shrink=0.8, aspect=30)
#     cbar.set_label("$|result|$", fontsize=12)
#     cbar.ax.tick_params(labelsize=10)
#     cbar.ax.yaxis.set_major_locator(LogLocator(base=10.0, numticks=6))
#     cbar.ax.yaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1, numticks=100))
#     cbar.ax.yaxis.set_major_formatter(LogFormatter(labelOnlyBase=False))
    
#     print(f"\nStatistics for gg={gg_value}:")
#     print(f"  Total points: {len(subset)}")
#     print(f"  Quenching points: {len(valid_subset)}")
#     print(f"  Stable points: {len(nan_subset)}")
#     if len(valid_subset) > 0:
#         print(f"  Quenching range: [{valid_subset['result'].min():.6f}, {valid_subset['result'].max():.6f}]")
#         print(f"  xs range: [{valid_subset['xs'].min():.2f}, {valid_subset['xs'].max():.2f}]")
#         print(f"  theta range: [{valid_subset['theta'].min():.2f}, {valid_subset['theta'].max():.2f}]")
    
#     plt.tight_layout()
#     plt.savefig(f"{output_dir}/single_gg_{gg_value}.pdf", dpi=300, bbox_inches='tight')
#     plt.show()

if __name__ == "__main__":
    file_path = "/globalhome/tus210/HPC/quenchin_actor/cpp_implementation/quench_new_actor_dynamic_box_early_update-4662663.out"   
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
    
    print("\n" + "="*50)
    print("Creating combined 2x2 subplot:")
    plot_heatmaps(df)
    
    # print("\n" + "="*50)
    # print("Creating four separate figures:")
    # plot_all_separate_figures(df)
    
    # You can also create plots without showing the data points:
    # plot_heatmaps(df, show_points=False)
    # plot_all_separate_figures(df, show_points=False)
