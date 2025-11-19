import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata
import matplotlib.ticker as ticker
import matplotlib.gridspec as gridspec
import gc

plt.close('all')  # Close previous figures
gc.collect()  # Force garbage collection

# File paths (modify to include more files if needed)
file_paths = [
     "TENO_HLLC_grid_906.dat"
    # ,"WENO_HLLC_grid_806.dat"
    # ,"MUSC_HLLC_grid_806.dat"
]

# Function to read data without pandas
def read_file(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
        data = [list(map(float, line.split())) for line in lines]
    return zip(*data)  # Return transposed data for easier plotting

# Create figure with GridSpec for better layout control
fig = plt.figure(figsize=(10, 6))
gs = gridspec.GridSpec(3, 4, height_ratios=[0.04, 1, 1], width_ratios=[1, 1, 1, 0.05])  # Colorbar row is now thinner

fig.suptitle("2D Riemann Problem: Density Distribution", fontsize=16, y=0.98)

labels = ['(a)', '(b)', '(c)']  # Labels for subplots
contour = None  # To store the last contour for colorbar

for i, file_path in enumerate(file_paths):
    x, y, density, *_ = read_file(file_path)
    grid_size = file_path.split('_')[0].replace('.dat', '').strip()  # Extract grid size

    if grid_size == 'WENO':
        grid_size = 'WENO-JS'
    elif grid_size == "MUSC":
        grid_size = 'MUSCL-THINC'

    # Define grid for contouring
    xi = np.linspace(min(x), 0.8, 1000)
    yi = np.linspace(min(y), 0.8, 1000)
    Xi, Yi = np.meshgrid(xi, yi)

    ax = fig.add_subplot(gs[1, i])  # Place in correct subplot grid position

    # Interpolate density data onto the grid
    Di = griddata((x, y), density, (Xi, Yi), method='cubic')

    # Plot filled contour in grayscale
    contour = ax.contourf(Xi, Yi, Di, levels=1000, cmap='gray')

    # Add subplot label (a), (b), (c)
    ax.text(0.05, 1.02, labels[i], transform=ax.transAxes, fontsize=12, fontweight='bold')
    ax.set_title(grid_size)
    ax.set_xlabel("x")
    ax.set_ylabel("y")

# Create the horizontal colorbar in the top-right corner
cax = fig.add_subplot(gs[0, 2])  # Colorbar in top-right corner
cbar = plt.colorbar(contour, cax=cax, orientation='horizontal')
cbar.formatter = ticker.FuncFormatter(lambda rho_tick, _: f'{rho_tick:.2f}')  # Adjust decimal places to 2
cbar.update_ticks()
cbar.ax.set_ylabel(r"$\rho$", rotation=0, labelpad=15, fontsize=14, verticalalignment='center')

plt.tight_layout(rect=[0, 0, 1, 0.98])  # Leave space for the suptitle
plt.subplots_adjust(hspace=0.3, wspace=0.3)  # Adjust spacing between plots and colorbar
plt.savefig("Euler.png")
plt.show()
