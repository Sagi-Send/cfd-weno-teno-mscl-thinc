import matplotlib.pyplot as plt
import numpy as np

# File paths
file_paths = [
     "TENO_HLLC_grid_610_Re_ 200.00.dat",
    "WENO_HLLC_grid_610_Re_ 200.00.dat",
    "MUSC_HLLC_grid_604_Re_ 200.00.dat"
]

file_path_reference = "reference4.csv"


# Function to read data without pandas
def read_file(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

        # Try to determine the delimiter
        data = []
        for line in lines:
            if ',' in line:  # Handle comma-separated values
                values = list(map(float, line.strip().split(',')))
            else:  # Handle space-separated values
                values = list(map(float, line.strip().split()))
            data.append(values)

    return np.array(data).T  # Transpose to get columns



# Create the first figure with custom layout
fig = plt.figure(figsize=(12, 6))
fig.suptitle("2D Viscous Shock Tube: Bottom Wall Density Distribution at "+r"$t=1$", fontsize=16, y=0.95)

colors = ["green", "red", "orange"]
count_colors = 0

# Main plot (left side, spans both rows)
for file_path in file_paths:
    data = read_file(file_path)
    x, y, density = data[0], data[1], data[2]

    scheme = file_path.split('_')[0].replace('.dat', '').strip()  # Extract grid size
    if scheme == 'WENO':
        scheme = 'WENO-JS'
        linest = "--"
    elif scheme == "MUSC":
        scheme = 'MUSCL-THINC'
        linest = "-."
    else:
        linest = ":"

    # Extract values where y is approximately 0.02
    x_plot = []
    rho_plot = []

    # y_target = 0.0505
    y_target = 0.0079
    tolerance = 1e-3  # Adjusted for accuracy
    for xi, yi, rho_i in zip(x, y, density):
        if np.isclose(yi, y_target, atol=tolerance) and xi >= 0.3:
            x_plot.append(xi)
            rho_plot.append(rho_i)

    # Sort data by x to ensure proper plotting
    sorted_indices = np.argsort(x_plot)
    x_plot = np.array(x_plot)[sorted_indices]
    rho_plot = np.array(rho_plot)[sorted_indices]

    plt.plot(x_plot, rho_plot, label=scheme, linewidth=1, linestyle=linest, color=colors[count_colors])
    count_colors += 1

# Plot reference
data = read_file(file_path_reference)
x, density = data[0], data[1]
plt.plot(x, density, label="Reference", linewidth=1, linestyle="-", color="black")

ax = plt.gca()  # Get current Axes
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)


plt.xlabel("x")
plt.ylabel("Density")
plt.legend()
plt.xlim(0.4, 0.97)

# Adjust layout and save the figure
plt.tight_layout(rect=[0, 0, 1, 0.96])  # Leave space for the suptitle
plt.subplots_adjust(hspace=0.45, wspace=0.2)  # Reduce space between rows and columns
plt.savefig("bottom_wall.png")
plt.show()
