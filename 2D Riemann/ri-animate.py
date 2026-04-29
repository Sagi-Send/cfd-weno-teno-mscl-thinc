import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np


# Snapshot files are written by ri-six.f90, one binary density stack per scheme.
SCHEMES = [
    ("TENO", "TENO", "(a)"),
    ("WENO", "WENO-JS", "(b)"),
    ("MUSC", "MUSCL-THINC", "(c)"),
]


def read_meta(path):
    # Metadata is a simple key=value text file paired with each binary stack.
    meta = {}
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line or "=" not in line:
                continue
            key, value = line.split("=", 1)
            meta[key] = value
    return meta


def frame_paths(frames_dir, scheme, grid_total, nsnap):
    # Keep the Python reader in lockstep with the Fortran snapshot naming.
    stem = f"rho_{scheme}_HLLC_grid_{grid_total}_frames_{nsnap}"
    return frames_dir / f"{stem}.bin", frames_dir / f"{stem}.meta"


def load_scheme(frames_dir, scheme, grid_total, nsnap):
    # Load and reshape the raw Fortran density dump into frame, y, x order.
    data_path, meta_path = frame_paths(frames_dir, scheme, grid_total, nsnap)
    if not data_path.exists():
        raise FileNotFoundError(f"Missing snapshot data: {data_path}")
    if not meta_path.exists():
        raise FileNotFoundError(f"Missing snapshot metadata: {meta_path}")

    meta = read_meta(meta_path)
    nx_phys = int(meta["nx_phys"])
    ny_phys = int(meta["ny_phys"])
    frame_count = int(meta["nsnap"])
    expected = frame_count * nx_phys * ny_phys
    data = np.fromfile(data_path, dtype=np.float64)
    if data.size != expected:
        raise ValueError(f"{data_path} has {data.size} values, expected {expected}")
    return data.reshape((frame_count, ny_phys, nx_phys)), meta


def main():
    parser = argparse.ArgumentParser(description="Render 2D Riemann density snapshots as an MP4.")
    parser.add_argument("--frames-dir", type=Path)
    parser.add_argument("--frames-dir-outflow", type=Path)
    parser.add_argument("--frames-dir-reflection", type=Path)
    parser.add_argument("--grid-total", required=True, type=int)
    parser.add_argument("--nsnap", required=True, type=int)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--ffmpeg-path", type=Path)
    parser.add_argument("--tfinal", type=float)
    parser.add_argument("--fps", type=int, default=12)
    parser.add_argument("--bitrate", type=int, default=2500)
    parser.add_argument("--crop-max", type=float, default=0.8)
    args = parser.parse_args()

    if args.ffmpeg_path:
        matplotlib.rcParams["animation.ffmpeg_path"] = str(args.ffmpeg_path)

    compare_mode = args.frames_dir_outflow is not None or args.frames_dir_reflection is not None
    if compare_mode:
        # Comparison mode stacks the same schemes for two boundary treatments.
        if args.frames_dir_outflow is None or args.frames_dir_reflection is None:
            raise ValueError("Both --frames-dir-outflow and --frames-dir-reflection are required for comparison plots")
        row_inputs = [
            ("Outflow", args.frames_dir_outflow),
            ("Reflection", args.frames_dir_reflection),
        ]
    else:
        if args.frames_dir is None:
            raise ValueError("--frames-dir is required unless comparison frame directories are provided")
        row_inputs = [("", args.frames_dir)]

    loaded_rows = []
    first_meta = None
    for row_label, frames_dir in row_inputs:
        loaded = []
        for scheme, title, label in SCHEMES:
            frames, meta = load_scheme(frames_dir, scheme, args.grid_total, args.nsnap)
            loaded.append((scheme, title, label, frames, meta))
            if first_meta is None:
                first_meta = meta
        loaded_rows.append((row_label, loaded))

    nx_phys = int(first_meta["nx_phys"])
    ny_phys = int(first_meta["ny_phys"])
    xl = float(first_meta["xl"])
    xr = float(first_meta["xr"])
    yl = float(first_meta["yl"])
    yr = float(first_meta["yr"])
    tfinal = args.tfinal if args.tfinal is not None else float(first_meta["tfinal"])

    x = np.linspace(xl + 0.5 * (xr - xl) / nx_phys, xr - 0.5 * (xr - xl) / nx_phys, nx_phys)
    y = np.linspace(yl + 0.5 * (yr - yl) / ny_phys, yr - 0.5 * (yr - yl) / ny_phys, ny_phys)
    x_mask = x <= args.crop_max
    y_mask = y <= args.crop_max
    extent = [x[x_mask][0], x[x_mask][-1], y[y_mask][0], y[y_mask][-1]]

    cropped_rows = []
    for row_label, loaded in loaded_rows:
        # Crop to the region used in the report figures.
        cropped = [(scheme, title, label, frames[:, y_mask, :][:, :, x_mask]) for scheme, title, label, frames, _ in loaded]
        cropped_rows.append((row_label, cropped))
    vmin = min(float(np.nanmin(frames)) for _, cropped in cropped_rows for _, _, _, frames in cropped)
    vmax = max(float(np.nanmax(frames)) for _, cropped in cropped_rows for _, _, _, frames in cropped)

    row_count = len(cropped_rows)
    fig, axes = plt.subplots(row_count, 3, figsize=(11.5, 4.1 * row_count), constrained_layout=False)
    if row_count == 1:
        axes = np.asarray([axes])
    fig.suptitle("2D Riemann Problem: Density Development", fontsize=15, y=0.96)
    images = []

    for row_idx, (row_label, cropped) in enumerate(cropped_rows):
        for col_idx, (ax, (_, title, label, frames)) in enumerate(zip(axes[row_idx], cropped)):
            # All images share one color scale so the schemes can be compared directly.
            image = ax.imshow(
                frames[0],
                origin="lower",
                extent=extent,
                cmap="gray",
                vmin=vmin,
                vmax=vmax,
                interpolation="bilinear",
                animated=True,
            )
            images.append((image, frames))
            if row_idx == 0:
                ax.set_title(title)
            ax.set_xlabel("x")
            ax.set_ylabel("y")
            if row_label and col_idx == 0:
                ax.text(-0.24, 0.5, row_label, transform=ax.transAxes, rotation=90, ha="center", va="center", fontsize=12)
            ax.text(0.04, 1.03, label, transform=ax.transAxes, fontsize=11, fontweight="bold")
            ax.set_xlim(extent[0], extent[1])
            ax.set_ylim(extent[2], extent[3])

    cbar = fig.colorbar(images[-1][0], ax=axes.ravel().tolist(), orientation="horizontal", fraction=0.05, pad=0.10)
    cbar.set_label(r"$\rho$", rotation=0)
    time_text = fig.text(0.5, 0.025, "", ha="center", va="center")

    def update(frame_idx):
        for image, frames in images:
            image.set_data(frames[frame_idx])
        time = tfinal * frame_idx / max(args.nsnap - 1, 1)
        time_text.set_text(f"t = {time:.4f}")
        return [*(image for image, _ in images), time_text]

    ani = animation.FuncAnimation(fig, update, frames=args.nsnap, interval=1000 / args.fps, blit=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    writer = animation.FFMpegWriter(fps=args.fps, bitrate=args.bitrate)
    ani.save(args.output, writer=writer, dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    main()
