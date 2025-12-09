#!/usr/bin/env python3
"""
plot_snapshots.py — visualize u1(x,t) snapshots produced by make_snapshots

Usage examples:
Single run visualization:
  python plot_snapshots.py --indir ./snap_out --preview 0
  python plot_snapshots.py --indir ./snap_out --png_dir frames --stride 2
  python plot_snapshots.py --indir ./snap_out --mp4 fhn.mp4 --fps 25
  python plot_snapshots.py --indir ./snap_out --panel grid.png --panel_rows 2 --panel_cols 3

Combined visualization (multiple A values stacked vertically):
  python plot_snapshots.py --indirs ./snap_A0_out ./snap_A052_out ./snap_A057_out --combined_png_dir combined_frames
  python plot_snapshots.py --indirs ./snap_A0_out ./snap_A052_out ./snap_A057_out --combined_mp4 combined.mp4

Notes:
- Uses matplotlib only (no seaborn). Requires ffmpeg available on PATH for MP4.
- Keeps y-limits fixed across frames for a stable animation.
- Combined mode stacks different A values vertically for easy comparison.
"""
import argparse, glob, json, os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter


def load_run(indir: str):
    meta_path = os.path.join(indir, 'meta.json')
    x_path = os.path.join(indir, 'x.csv')
    snap_files = sorted(glob.glob(os.path.join(indir, 'snap_*.csv')))
    if not os.path.exists(meta_path) or not os.path.exists(x_path) or not snap_files:
        raise FileNotFoundError("indir must contain meta.json, x.csv, and snap_*.csv")

    with open(meta_path) as f:
        meta = json.load(f)
    x = np.loadtxt(x_path, delimiter=',')
    U = np.stack([np.loadtxt(f, delimiter=',') for f in snap_files], axis=0)  # (K, N)
    return x, U, meta


def compute_limits(x, U, ypad=0.05):
    xmin, xmax = float(np.min(x)), float(np.max(x))
    umin, umax = float(np.min(U)), float(np.max(U))
    if umax == umin:
        umax = umin + 1.0
    pad = ypad * (umax - umin)
    return (xmin, xmax), (umin - pad, umax + pad)


def save_frames(indir: str, outdir: str, stride: int = 1, dpi: int = 120,
                title_extra: str = ""):
    os.makedirs(outdir, exist_ok=True)
    x, U, meta = load_run(indir)
    K, N = U.shape
    (xmin, xmax), (ymin, ymax) = compute_limits(x, U)

    t0, tf = meta['t0'], meta['tf']
    for k in range(0, K, max(1, stride)):
        t = t0 + (tf - t0) * (k / (K - 1)) if K > 1 else t0
        fig = plt.figure(figsize=(8, 3))
        ax = fig.add_subplot(111)
        ax.plot(x, U[k])
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        ax.set_xlabel('x')
        ax.set_ylabel('u(x,t)')
        base = f"A={meta['A']}"
        ax.set_title(f"{base}  |  t={t:.2f} {title_extra}")
        fig.tight_layout()
        fig.savefig(os.path.join(outdir, f"frame_{k:04d}.png"), dpi=dpi)
        plt.close(fig)


def save_mp4(indir: str, mp4_path: str, fps: int = 25, dpi: int = 120,
             title_extra: str = ""):
    x, U, meta = load_run(indir)
    K, N = U.shape
    (xmin, xmax), (ymin, ymax) = compute_limits(x, U)

    fig = plt.figure(figsize=(8, 3))
    ax = fig.add_subplot(111)
    (line,) = ax.plot([], [])
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_xlabel('x')
    ax.set_ylabel('u(x,t)')

    base = f"A={meta['A']}"
    t0, tf = meta['t0'], meta['tf']

    def update(k):
        t = t0 + (tf - t0) * (k / (K - 1)) if K > 1 else t0
        line.set_data(x, U[k])
        ax.set_title(f"FHN wave  |  {base}  |  t={t:.2f} {title_extra}")
        return (line,)

    writer = FFMpegWriter(fps=fps, metadata=dict(artist='FHN'))
    with writer.saving(fig, mp4_path, dpi=dpi):
        for k in range(K):
            update(k)
            writer.grab_frame()
    plt.close(fig)


def save_panel(indir: str, out_png: str, rows: int, cols: int, start: int = 0, stride: int = 1, dpi: int = 150):
    """Create a rows×cols panel of snapshots for slides."""
    x, U, meta = load_run(indir)
    K, N = U.shape
    (xmin, xmax), (ymin, ymax) = compute_limits(x, U)
    t0, tf = meta['t0'], meta['tf']

    total = rows * cols
    idxs = [start + i * stride for i in range(total)]
    idxs = [min(k, K - 1) for k in idxs]

    fig = plt.figure(figsize=(cols * 4.5, rows * 3.2))
    for i, k in enumerate(idxs):
        ax = fig.add_subplot(rows, cols, i + 1)
        t = t0 + (tf - t0) * (k / (K - 1)) if K > 1 else t0
        ax.plot(x, U[k])
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        ax.set_title(f"t={t:.2f}")
        if i // cols == rows - 1:
            ax.set_xlabel('x')
        if i % cols == 0:
            ax.set_ylabel('u₁(x,t)')
    fig.suptitle(f"FHN wave  •  gg={meta['gg']}, θ={meta['theta']}, xₛ={meta['xs']}, A={meta['A']}")
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    fig.savefig(out_png, dpi=dpi)
    plt.close(fig)


def save_combined_frames(indirs: list, outdir: str, stride: int = 1, dpi: int = 120):
    """Create combined frames stacking multiple simulation results vertically."""
    os.makedirs(outdir, exist_ok=True)
    
    # Load all runs
    runs = []
    for indir in indirs:
        x, U, meta = load_run(indir)
        runs.append((x, U, meta))
    
    # Find common limits across all runs
    all_x = np.concatenate([x for x, _, _ in runs])
    all_U = np.concatenate([U for _, U, _ in runs])
    (xmin, xmax), (ymin, ymax) = compute_limits(all_x, all_U)
    
    # Use the minimum number of frames across all runs
    min_frames = min(U.shape[0] for _, U, _ in runs)
    
    for k in range(0, min_frames, max(1, stride)):
        fig, axes = plt.subplots(len(runs), 1, figsize=(10, 3 * len(runs)), sharex=True)
        if len(runs) == 1:
            axes = [axes]  # Make it a list for consistency
        
        for i, (x, U, meta) in enumerate(runs):
            t0, tf = meta['t0'], meta['tf']
            t = t0 + (tf - t0) * (k / (U.shape[0] - 1)) if U.shape[0] > 1 else t0
            
            axes[i].plot(x, U[k])
            axes[i].set_xlim(xmin, xmax)
            axes[i].set_ylim(ymin, ymax)
            axes[i].set_ylabel('u(x,t)')
            axes[i].set_title(f"A={meta['A']}, t={t:.2f}")
            axes[i].grid(True, alpha=0.3)
        
        # Only add x-label to bottom subplot
        axes[-1].set_xlabel('x')
        
        fig.tight_layout()
        fig.savefig(os.path.join(outdir, f"combined_frame_{k:04d}.png"), dpi=dpi)
        plt.close(fig)


def save_combined_mp4(indirs: list, mp4_path: str, fps: int = 25, dpi: int = 120):
    """Create MP4 animation with multiple simulations stacked vertically."""
    # Load all runs
    runs = []
    for indir in indirs:
        x, U, meta = load_run(indir)
        runs.append((x, U, meta))
    
    # Find common limits across all runs
    all_x = np.concatenate([x for x, _, _ in runs])
    all_U = np.concatenate([U for _, U, _ in runs])
    (xmin, xmax), (ymin, ymax) = compute_limits(all_x, all_U)
    
    # Use the minimum number of frames across all runs
    min_frames = min(U.shape[0] for _, U, _ in runs)
    
    fig, axes = plt.subplots(len(runs), 1, figsize=(10, 3 * len(runs)), sharex=True)
    if len(runs) == 1:
        axes = [axes]  # Make it a list for consistency
    
    lines = []
    for i, (x, U, meta) in enumerate(runs):
        line, = axes[i].plot([], [])
        lines.append(line)
        axes[i].set_xlim(xmin, xmax)
        axes[i].set_ylim(ymin, ymax)
        axes[i].set_ylabel('u(x,t)')
        axes[i].grid(True, alpha=0.3)
    
    axes[-1].set_xlabel('x')
    
    def update(k):
        for i, (x, U, meta) in enumerate(runs):
            t0, tf = meta['t0'], meta['tf']
            t = t0 + (tf - t0) * (k / (U.shape[0] - 1)) if U.shape[0] > 1 else t0
            lines[i].set_data(x, U[k])
            axes[i].set_title(f"A={meta['A']}, t={t:.2f}")
        return lines

    writer = FFMpegWriter(fps=fps, metadata=dict(artist='FHN'))
    with writer.saving(fig, mp4_path, dpi=dpi):
        for k in range(min_frames):
            update(k)
            writer.grab_frame()
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--indir', help='folder with meta.json, x.csv, snap_*.csv (for single run)')
    ap.add_argument('--indirs', nargs='+', help='multiple folders for combined visualization')
    ap.add_argument('--preview', type=int, default=-1, help='show a single frame index and exit')
    ap.add_argument('--png_dir', default='', help='if set, save individual frames here')
    ap.add_argument('--stride', type=int, default=1, help='frame stride for --png_dir')
    ap.add_argument('--mp4', default='', help='if set, write MP4 animation to this path')
    ap.add_argument('--fps', type=int, default=25, help='frames per second for MP4')
    ap.add_argument('--panel', default='', help='if set, save a rows×cols grid panel to this PNG')
    ap.add_argument('--panel_rows', type=int, default=2)
    ap.add_argument('--panel_cols', type=int, default=3)
    ap.add_argument('--panel_start', type=int, default=0)
    ap.add_argument('--panel_stride', type=int, default=1)
    ap.add_argument('--combined_png_dir', default='', help='if set with --indirs, save combined frames here')
    ap.add_argument('--combined_mp4', default='', help='if set with --indirs, save combined MP4 here')
    args = ap.parse_args()
    
    # Check if we have directories specified
    if not args.indir and not args.indirs:
        ap.error("Either --indir (single run) or --indirs (multiple runs) must be specified")
    
    # Handle combined visualization (multiple directories)
    if args.indirs:
        if args.combined_png_dir:
            save_combined_frames(args.indirs, args.combined_png_dir, stride=max(1, args.stride))
            print(f"Combined frames saved to {args.combined_png_dir}")
        
        if args.combined_mp4:
            save_combined_mp4(args.indirs, args.combined_mp4, fps=args.fps)
            print(f"Combined MP4 saved to {args.combined_mp4}")
        
        # No other operations supported with --indirs
        return
    
    # Single directory operations (existing functionality)
    if not args.indir:
        ap.error("--indir must be specified for single run operations")

    # preview single frame
    if args.preview >= 0:
        x, U, meta = load_run(args.indir)
        (xmin, xmax), (ymin, ymax) = compute_limits(x, U)
        k = min(max(0, args.preview), U.shape[0]-1)
        t0, tf = meta['t0'], meta['tf']
        t = t0 + (tf - t0) * (k / (U.shape[0] - 1)) if U.shape[0] > 1 else t0
        plt.figure(figsize=(8,3))
        plt.plot(x, U[k])
        plt.xlim(xmin, xmax); plt.ylim(ymin, ymax)
        plt.xlabel('x'); plt.ylabel('u₁(x,t)')
        base = f"gg={meta['gg']}, θ={meta['theta']}, xₛ={meta['xs']}, A={meta['A']}"
        plt.title(f"FHN wave  |  {base}  |  t={t:.2f}")
        plt.tight_layout(); plt.show()
        return

    if args.png_dir:
        save_frames(args.indir, args.png_dir, stride=max(1, args.stride))

    if args.mp4:
        save_mp4(args.indir, args.mp4, fps=args.fps)

    if args.panel:
        save_panel(args.indir, args.panel, rows=args.panel_rows, cols=args.panel_cols,
                   start=args.panel_start, stride=max(1, args.panel_stride))


if __name__ == '__main__':
    main()
