import matplotlib.pyplot as plt
import numpy as np

# File paths
file_paths = [
     "TENO.dat"
]


# Function to read data without pandas
def read_file(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
        data = [list(map(float, line.split())) for line in lines]
    return zip(*data)  # Return transposed data for easier plotting


def exact_solution(x_exact, t=2.0, u=0.1):
    return 1 + 0.98 * np.sin(2 * np.pi * (x_exact - u * t))


def initial_condition(x_initial, t=0, u=0.1):
    return 1 + 0.98 * np.sin(2 * np.pi * (x_exact - u * t))


# Create the figure with custom layout
fig = plt.figure(figsize=(12, 6))
fig.suptitle("Density Wave Test Case", fontsize=16, y=0.95)

# Add exact solution
x_exact = np.linspace(-1, 1, 500)
y_exact = exact_solution(x_exact)


# -------------------------------------
# 2. Plot the L1 norm vs. grid spacing
# -------------------------------------

# Read the L1 norm data
grid_sizes, l1 = read_file("TENO.dat")
grid_sizes = np.array(grid_sizes)
l1 = np.array(l1)

# Convert grid sizes to h = 1 / grid_sizes (or whichever measure is appropriate)
# If "grid_sizes" in your data file already represents number of cells, you can do:
# h = 1 / grid_sizes
# or, if "grid_sizes" is already the grid spacing, use them directly.
# Here, we'll assume grid_sizes is the number of cells, and h = 1 / grid_sizes:
h = 1.0 / grid_sizes

ax_sec = plt.subplot2grid((1, 2), (0, 1))

ax_sec.plot(h, l1, marker='o', label='TENO', linestyle='-', color='orange')

# grid_sizes, l1 = read_file("WENO.dat")
grid_sizes = np.array(grid_sizes)
l1 = np.array(l1)

# Convert grid sizes to h = 1 / grid_sizes (or whichever measure is appropriate)
# If "grid_sizes" in your data file already represents number of cells, you can do:
# h = 1 / grid_sizes
# or, if "grid_sizes" is already the grid spacing, use them directly.
# Here, we'll assume grid_sizes is the number of cells, and h = 1 / grid_sizes:
h = 1.0 / grid_sizes


ax_sec.plot(h, l1, marker='o', label='WENO', linestyle='-', color='orange')

# -----------------------------------------------------------
# 3. Add a second-order reference line (slope ~ 2 in log-log)
# -----------------------------------------------------------
# We'll reference the last data point so the line matches the scale of our plot.
# If you want it to pass through the first point, you can adjust accordingly.

h_ref = h[-1]      # choose the last h
l1_ref = l1[-1]    # corresponding L1 norm
slope = 2.0        # second-order slope

# In log-log space: log(L1) = log(C) + slope*log(h)
# => L1 = C * h^slope
# We solve for C: C = l1_ref / (h_ref^slope)
C = l1_ref / (h_ref**slope)

# Generate a smooth array of h-values for the reference line
h_line = np.logspace(np.log10(h[1]), np.log10(h[-1]), 10)
l1_line = C * (h_line**slope)

ax_sec.plot(h_line, l1_line, '--', label='2nd Order Reference')

# Customize, scale, and save
ax_sec.grid(False)
ax_sec.set_xscale('log')
ax_sec.set_yscale('log')
ax_sec.set_xlabel("Grid Spacing (h)")
ax_sec.set_ylabel("L1 Norm of Error")
ax_sec.set_title("L1 Norm of Error as a Function of Grid Resolution")
ax_sec.legend()

# Adjust layout and save the figure
plt.tight_layout(rect=[0, 0, 1, 0.96])  # Leave space for the title
plt.savefig("51.png")
plt.show()
