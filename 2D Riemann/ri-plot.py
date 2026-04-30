from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MaxNLocator


# Plot the 2D Riemann density fields used in the README.
FRAME_DIR = Path("frames_t1p5")
GRID_TOTAL = 806
SNAPSHOTS = 100
TARGET_TIME = 0.8

CASES = [
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
            if line and "=" in line:
                key, value = line.split("=", 1)
                meta[key] = value.strip()
    return meta


def snapshot_paths(scheme):
    # Keep the figure generator in lockstep with the Fortran snapshot naming.
    stem = f"rho_{scheme}_HLLC_grid_{GRID_TOTAL}_frames_{SNAPSHOTS}"
    return FRAME_DIR / f"{stem}.bin", FRAME_DIR / f"{stem}.meta"


def read_density_frame(scheme):
    # Read only the frame nearest TARGET_TIME from the y-major binary stack.
    data_path, meta_path = snapshot_paths(scheme)
    meta = read_meta(meta_path)
    nx_phys = int(meta["nx_phys"])
    ny_phys = int(meta["ny_phys"])
    nsnap = int(meta["nsnap"])
    final_time = float(meta["tfinal"])
    frame_index = int(round(TARGET_TIME / final_time * (nsnap - 1)))
    frame_index = max(0, min(frame_index, nsnap - 1))
    frame_time = final_time * frame_index / max(nsnap - 1, 1)
    frame_size = nx_phys * ny_phys
    offset = frame_index * frame_size * np.dtype(np.float64).itemsize
    frame = np.memmap(data_path, dtype=np.float64, mode="r", offset=offset, shape=(ny_phys, nx_phys))

    xl = float(meta["xl"])
    xr = float(meta["xr"])
    yl = float(meta["yl"])
    yr = float(meta["yr"])
    x = np.linspace(xl + 0.5 * (xr - xl) / nx_phys, xr - 0.5 * (xr - xl) / nx_phys, nx_phys)
    y = np.linspace(yl + 0.5 * (yr - yl) / ny_phys, yr - 0.5 * (yr - yl) / ny_phys, ny_phys)
    return x, y, np.asarray(frame), frame_time


def crop_to_report_window(x, y, density, crop_max=0.8):
    # Match the report view and keep the README image focused on the interaction zone.
    x_mask = x <= crop_max
    y_mask = y <= crop_max
    cropped = density[np.ix_(y_mask, x_mask)]
    extent = [x[x_mask][0], x[x_mask][-1], y[y_mask][0], y[y_mask][-1]]
    return cropped, extent


def main():
    loaded = []
    final_time = None
    for scheme, title, label in CASES:
        x, y, density, final_time = read_density_frame(scheme)
        cropped, extent = crop_to_report_window(x, y, density)
        loaded.append((title, label, cropped, extent))

    vmin = min(float(np.nanmin(density)) for _, _, density, _ in loaded)
    vmax = max(float(np.nanmax(density)) for _, _, density, _ in loaded)

    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.6), constrained_layout=True)
    fig.suptitle(f"2D Riemann Problem: Density Distribution at t = {TARGET_TIME:.1f}", fontsize=18)

    image = None
    for ax, (title, label, density, extent) in zip(axes, loaded):
        image = ax.imshow(
            density,
            origin="lower",
            extent=extent,
            cmap="gray",
            vmin=vmin,
            vmax=vmax,
            interpolation="bilinear",
            aspect="equal",
        )
        ax.set_title(title, fontsize=14)
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.text(0.04, 1.03, label, transform=ax.transAxes, fontsize=12, fontweight="bold")
        ax.set_xlim(extent[0], extent[1])
        ax.set_ylim(extent[2], extent[3])

    # Put the shared legend below the panels with sparse ticks to avoid label collisions.
    colorbar = fig.colorbar(image, ax=axes, orientation="horizontal", fraction=0.08, pad=0.14)
    colorbar.locator = MaxNLocator(nbins=6)
    colorbar.update_ticks()
    colorbar.ax.tick_params(labelsize=10)
    colorbar.set_label(r"$\rho$", rotation=0, fontsize=13)

    fig.savefig("Euler.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
