import matplotlib.pyplot as plt
import numpy as np

# Load the data from the file
file_path = 'SodShockTube.dat'
data = np.loadtxt(file_path, skiprows=1)

# Extract the required columns
x_exact = data[:, 0] + 0.5  # First column

density_exact = data[:, 1]  # Second column
velocity_exact = data[:, 2]  # Third column
pressure_exact = data[:, 3]  # Fourth column

# File paths
file_paths = [
     "TENO_HLLC_grid_806.dat"
    ,"MUSC_HLLC_grid_806.dat"
    ,"WENO_HLLC_grid_806.dat"
]


# Function to read data without pandas
def read_file(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
        data = [list(map(float, line.split())) for line in lines]
    return zip(*data)  # Return transposed data for easier plotting


# Create the first figure with custom layout
fig = plt.figure(figsize=(12, 6))
fig.suptitle("Sod's Shock Tube Problem: Density Distribution", fontsize=16, y=0.95)

# Main plot (left side, spans both rows)
ax_main = plt.subplot2grid((2, 2), (0, 0), rowspan=2)
for file_path in file_paths:
    x, y, *_ = read_file(file_path)
    grid_size = file_path.split('_')[0].replace('.dat', '').strip()  # Extract grid size
    if grid_size == 'WENO':
        grid_size = 'WENO-JS'
        linest = "--"
    elif grid_size == "MUSC":
        grid_size = 'MUSCL-THINC'
        linest = "-."
    else:
        linest = "-"

    grid_size = grid_size
    ax_main.plot(x, y, label=grid_size, linewidth=1, linestyle=linest)
ax_main.plot(x_exact, density_exact, label='Exact', linestyle=':', color='blue')
ax_main.set_title("Overall Density Distribution")
ax_main.set_xlabel("x")
ax_main.set_ylabel("Density")
ax_main.legend()
ax_main.set_xlim(0, 1)

# Zoomed-in plot 1 (top-right)
ax_zoom1 = plt.subplot2grid((2, 2), (0, 1))
for file_path in file_paths:
    x, y, *_ = read_file(file_path)
    grid_size = file_path.split('_')[0].replace('.dat', '').strip()  # Extract grid size
    solver = file_path.split('_')[1].replace('.dat', '').strip()
    if grid_size == 'WENO':
        grid_size = 'WENO-JS'
        linest="--"
    elif grid_size == "MUSC":
        grid_size = 'MUSCL-THINC'
        linest = "-."
    else:
        linest = "-"
    grid_size = grid_size
    ax_zoom1.plot(x, y, label=grid_size, linewidth=1, linestyle=linest)
ax_zoom1.plot(x_exact, density_exact, label='Exact', linestyle=':', color='blue')
ax_zoom1.set_title("Zoomed-in View: Density Discontinuity")
ax_zoom1.set_xlabel("x")
ax_zoom1.set_ylabel("Density")
ax_zoom1.set_xlim(0.67, 0.7)
ax_zoom1.set_ylim(0.25, 0.45)
ax_zoom1.legend()

# Zoomed-in plot 2 (bottom-right)
ax_zoom2 = plt.subplot2grid((2, 2), (1, 1))
for file_path in file_paths:
    x, y, *_ = read_file(file_path)
    grid_size = file_path.split('_')[0].replace('.dat', '').strip()  # Extract grid size
    solver = file_path.split('_')[1].replace('.dat', '').strip()
    if grid_size == 'WENO':
        grid_size = 'WENO-JS'
        linest="--"
    elif grid_size == "MUSC":
        grid_size = 'MUSCL-THINC'
        linest = "-."
    else:
        linest = "-"
    grid_size = grid_size
    ax_zoom2.plot(x, y, label=grid_size, linewidth=1, linestyle=linest)
ax_zoom2.plot(x_exact, density_exact, label='Exact', linestyle=':', color='blue')
ax_zoom2.set_title("Zoomed-in View: Shock Wave")
ax_zoom2.set_xlabel("x")
ax_zoom2.set_ylabel("Density")
ax_zoom2.set_xlim(0.846, 0.855)
ax_zoom2.set_ylim(0.1, 0.275)
ax_zoom2.legend()

# Adjust layout and save the figure
plt.tight_layout(rect=[0, 0, 1, 0.96])  # Leave space for the suptitle
plt.subplots_adjust(hspace=0.45, wspace=0.15)  # Reduce space between rows and columns
plt.savefig("Sod Tube")
plt.show()
