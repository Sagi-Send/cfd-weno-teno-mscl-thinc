import re
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata
import multiprocessing

# File paths: modify this array as needed.
file_paths = [
    "TENO_HLLC_grid_610_Re_ 200.00.dat",
    "WENO_HLLC_grid_610_Re_ 200.00.dat",
    "MUSC_HLLC_grid_604_Re_ 200.00.dat",

    "TENO_HLLC_grid_610_Re_1000.00.dat",
    "WENO_HLLC_grid_610_Re_1000.00.dat",
    "MUSC_HLLC_grid_604_Re_1000.00.dat"
]

# Grid resolution and domain
grid_res = 1000
xi = np.linspace(0.35, 0.98, grid_res)
yi = np.linspace(0.02, 0.25, grid_res)
Xi, Yi = np.meshgrid(xi, yi)

# y threshold for switching methods
y_threshold = 0

# --- Group files by Reynolds number ---
groups = {}
for fp in file_paths:
    m = re.search(r"Re_\s*([\d\.]+)", fp)
    if m:
        re_val = m.group(1).rstrip('.')  # Remove trailing period if present.
    else:
        re_val = "unknown"
    groups.setdefault(re_val, []).append(fp)

# Sort groups numerically (if possible)
sorted_re_keys = sorted(groups.keys(), key=lambda x: float(x) if x != "unknown" else float('inf'))

# Grid dimensions: one row per unique Re; columns = max files per group
nrows = len(sorted_re_keys) + 1
ncols = max(len(groups[re_val]) for re_val in sorted_re_keys)

# Build tasks for multiprocessing
tasks = []
for row_idx, re_val in enumerate(sorted_re_keys):
    group_files = groups[re_val]
    for col_idx, fp in enumerate(group_files):
        tasks.append((row_idx, col_idx, fp))

def read_file(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
        data = [list(map(float, line.split())) for line in lines]
    return np.array(data).T

def piecewise_interpolation(x, y, values, Xi, Yi, threshold=0.03):
    """
    Interpolate linearly where Y < threshold, and cubically otherwise.
    Return the combined array.
    """
    # 1) Interpolate with linear
    linear_vals = griddata((x, y), values, (Xi, Yi), method='linear')
    # 2) Interpolate with cubic
    cubic_vals = griddata((x, y), values, (Xi, Yi), method='cubic')

    # Create a combined result array
    Di = np.zeros_like(linear_vals)

    # Define masks
    mask_linear = (Yi < threshold)
    mask_cubic  = ~mask_linear  # everything else

    # Combine: use linear interpolation where y < threshold, cubic elsewhere
    Di[mask_linear] = linear_vals[mask_linear]
    Di[mask_cubic]  = cubic_vals[mask_cubic]

    return Di

def process_subplot(args):
    row, col, fp = args
    data = read_file(fp)
    x, y, density, temp = data[0], data[1], data[2], data[3]

    # Interpolate density field then compute its gradient magnitude
    Di = piecewise_interpolation(x, y, density, Xi, Yi, threshold=y_threshold)
    grad_x = np.gradient(Di, xi, axis=1)
    grad_y = np.gradient(Di, yi, axis=0)
    grad_mag = np.hypot(grad_x, grad_y)
    return row, col, grad_mag, fp

def process_subplot_temp(args):
    row, col, fp = args
    data = read_file(fp)
    x, y, density, temp = data[0], data[1], data[2], data[3]

    # Interpolate the temperature field directly
    temp_field = piecewise_interpolation(x, y, temp, Xi, Yi, threshold=y_threshold)
    return row, col, temp_field, fp

if __name__ == '__main__':
    # ----------------- Density Gradient Plot -----------------
    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 4 * nrows))
    fig.suptitle("2D Viscous Shock Tube: Density Gradient Distribution at "+r"$t=1$",
                 fontsize=16, y=0.97)

    # Ensure axes is a 2D array
    if nrows == 1:
        axes = np.expand_dims(axes, axis=0)
    if ncols == 1:
        axes = np.expand_dims(axes, axis=1)

    # Process files in parallel for density gradient
    with multiprocessing.Pool() as pool:
        results = pool.map(process_subplot, tasks)

    # Plot each subplot for density gradient
    plotted = {}
    for row, col, grad_mag, fp in results:
        ax = axes[row, col]
        c = ax.contourf(Xi, Yi, grad_mag, levels=1000, cmap='gray_r')
        ax.set_xlabel("x")
        ax.set_ylabel("y")

        # Derive title from file naming scheme
        scheme = fp.split('_')[0]
        if scheme == "TENO-HLLC":
            scheme_title = "TENO"
        elif scheme == "WENO-HLLC" or scheme == "WENO":
            scheme_title = "WENO-JS"
        elif scheme == "MUSC":
            scheme_title = "MUSCL-THINC"
        else:
            scheme_title = scheme
        ax.set_title(scheme_title)
        plotted[(row, col)] = True

    # Turn off any unused subplots for density gradient
    for row in range(nrows):
        for col in range(ncols):
            if (row, col) not in plotted:
                axes[row, col].axis('off')

    # Label each row with (a), (b), etc.
    for row in range(nrows):
        label = f"({chr(97 + row)})"
        ax_left = axes[row, 0]
        ax_left.text(-0.1, 0.95, label,
                     transform=ax_left.transAxes,
                     ha='right', va='center',
                     fontsize=14, fontweight='bold')

    plt.tight_layout(rect=[0, 0, 1, 0.94])
    plt.subplots_adjust(left=0.05, right=0.99, hspace=0.4, wspace=0.15)

    # Save the density gradient figure (choose appropriate time label)
    plt.savefig("Gradient1.png", dpi=600)

    # ----------------- Temperature Plot -----------------
    fig_temp, axes_temp = plt.subplots(nrows, ncols, figsize=(4 * ncols, 4 * nrows))
    fig_temp.suptitle("2D Viscous Shock Tube: Temperature Distribution at "+r"$t=1$", fontsize=16, y=0.97)

    # Ensure axes_temp is 2D
    if nrows == 1:
        axes_temp = np.expand_dims(axes_temp, axis=0)
    if ncols == 1:
        axes_temp = np.expand_dims(axes_temp, axis=1)

    # Process files in parallel for temperature field
    with multiprocessing.Pool() as pool:
        temp_results = pool.map(process_subplot_temp, tasks)

    # Plot each subplot for temperature field
    plotted_temp = {}
    for row, col, temp_field, fp in temp_results:
        ax_temp = axes_temp[row, col]
        c_temp = ax_temp.contourf(Xi, Yi, temp_field, levels=1000, cmap='gray_r')
        ax_temp.set_xlabel("x")
        ax_temp.set_ylabel("y")

        # Derive title from file naming scheme
        scheme = fp.split('_')[0]
        if scheme == "TENO-HLLC":
            scheme_title = "TENO"
        elif scheme == "WENO-HLLC" or scheme == "WENO":
            scheme_title = "WENO-JS"
        elif scheme == "MUSC":
            scheme_title = "MUSCL-THINC"
        else:
            scheme_title = scheme
        ax_temp.set_title(scheme_title)
        plotted_temp[(row, col)] = True

    # Turn off any unused subplots for temperature plot
    for row in range(nrows):
        for col in range(ncols):
            if (row, col) not in plotted_temp:
                axes_temp[row, col].axis('off')

    # Label each row with (a), (b), etc.
    for row in range(nrows):
        label = f"({chr(97 + row)})"
        ax_left = axes_temp[row, 0]
        ax_left.text(-0.1, 0.95, label,
                     transform=ax_left.transAxes,
                     ha='right', va='center',
                     fontsize=14, fontweight='bold')

    plt.tight_layout(rect=[0, 0, 1, 0.94])
    plt.subplots_adjust(left=0.05, right=0.99, hspace=0.4, wspace=0.15)

    # Save the temperature figure
    plt.savefig("Temperature.png", dpi=600)

    # Show both figures
    plt.show()
