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
     "TENO_HLLC_ct_-2.dat"
    ,"TENO_HLLC_ct_-4.dat"
    ,"TENO_HLLC_ct_-5.dat"
    ,"TENO_HLLC_ct_-8.dat"
    ,"MUSC_HLLC_beta_1.6.dat"
    ,"MUSC_HLLC_beta_2.0.dat"
    ,"MUSC_HLLC_beta_2.4.dat"
    ,"MUSC_HLLC_beta_2.8.dat"
]

linestyles=[":", "--", "-.", "-"]
# Function to read data without pandas
def read_file(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
        data = [list(map(float, line.split())) for line in lines]
    return zip(*data)  # Return transposed data for easier plotting


# Create the first figure with custom layout
fig = plt.figure(figsize=(12, 6))
fig.suptitle("MUSCL-THINC " + r"$\beta$" + " Paramteric Study: Sod's Shock Tube Density Distribution", fontsize=16, y=0.95)


# MUSCL-THINC shock (top-left)
ax_zoom1 = plt.subplot2grid((2, 2), (0, 0))
markers = ['o', 's', '<', '^']
count = 0
for file_path in file_paths:
    x, y, *_ = read_file(file_path)
    scheme = file_path.split('_')[0].replace('.dat', '').strip()  # Extract grid size
    if scheme == "MUSC":
        param = file_path.split('_')[3].replace('.dat', '').strip()  # Extract grid size
        ax_zoom1.plot(x, y, label=r"$\beta$="+param, linewidth=0.8, marker=markers[count], markersize=4, linestyle=linestyles[count])
        count += 1
ax_zoom1.plot(x_exact, density_exact, label='Exact', linestyle=':', color='blue')
ax_zoom1.set_title("Density Discontinuity")
ax_zoom1.set_xlabel("x")
ax_zoom1.set_ylabel("Density")
ax_zoom1.set_xlim(0.675, 0.695)
ax_zoom1.set_ylim(0.25, 0.45)
ax_zoom1.legend()
ax_zoom1.spines['top'].set_visible(False)
ax_zoom1.spines['right'].set_visible(False)

# MUSCL-THINC density discontinuity (bottom-left)
ax_zoom2 = plt.subplot2grid((2, 2), (0, 1))
count = 0
for file_path in file_paths:
    x, y, *_ = read_file(file_path)
    scheme = file_path.split('_')[0].replace('.dat', '').strip()  # Extract grid size
    if scheme == "MUSC":
        param = file_path.split('_')[3].replace('.dat', '').strip()  # Extract grid size
        ax_zoom2.plot(x, y, label=r"$\beta$="+param, linewidth=0.8, marker=markers[count], markersize=4, linestyle=linestyles[count])
        count += 1
ax_zoom2.plot(x_exact, density_exact, label='Exact', linestyle=':', color='blue')
ax_zoom2.set_title("Shock Wave")
ax_zoom2.set_xlabel("x")
ax_zoom2.set_ylabel("Density")
ax_zoom2.set_xlim(0.846, 0.855)
ax_zoom2.set_ylim(0.1, 0.275)
ax_zoom2.legend()
ax_zoom2.spines['top'].set_visible(False)
ax_zoom2.spines['right'].set_visible(False)

# Adjust layout and save the figure
plt.tight_layout(rect=[0, 0, 1, 0.96])  # Leave space for the suptitle
plt.subplots_adjust(hspace=0.45, wspace=0.15)  # Reduce space between rows and columns
plt.savefig("Calibration-MUSCL")
plt.show()



# Create the first figure with custom layout
fig = plt.figure(figsize=(12, 6))
fig.suptitle("TENO " + r"$C_T$" + " Paramteric Study: Sod's Shock Tube Density Distribution", fontsize=16, y=0.95)

# TENO density discontinuity (top-right)
ax_zoom3 = plt.subplot2grid((2, 2), (0, 0))
count = 0
for file_path in file_paths:
    x, y, *_ = read_file(file_path)
    scheme = file_path.split('_')[0].replace('.dat', '').strip()  # Extract grid size
    if scheme == "TENO":
        param = file_path.split('_')[3].replace('.dat', '').strip()  # Extract grid size
        if param == "-5":
            param = "-7"
        ax_zoom3.plot(x, y, label=r"$C_T = 10^{" + str(param) + "}$", linewidth=0.8, marker=markers[count], markersize=4, linestyle=linestyles[count])
        count += 1
ax_zoom3.plot(x_exact, density_exact, label='Exact', linestyle=':', color='blue')
ax_zoom3.set_title("Density Discontinuity")
ax_zoom3.set_xlabel("x")
ax_zoom3.set_ylabel("Density")
ax_zoom3.set_xlim(0.675, 0.695)
ax_zoom3.set_ylim(0.25, 0.45)
ax_zoom3.legend()
ax_zoom3.spines['top'].set_visible(False)
ax_zoom3.spines['right'].set_visible(False)


# TENO shock (bottom-right)
ax_zoom4 = plt.subplot2grid((2, 2), (0, 1))
count = 0
for file_path in file_paths:
    x, y, *_ = read_file(file_path)
    scheme = file_path.split('_')[0].replace('.dat', '').strip()  # Extract grid size
    if scheme == "TENO":
        param = file_path.split('_')[3].replace('.dat', '').strip()  # Extract grid size
        if param == "-5":
            param = "-7"
        ax_zoom4.plot(x, y, label=r"$C_T = 10^{" + str(param) + "}$", linewidth=0.8, marker=markers[count], markersize=4, linestyle=linestyles[count])
        count += 1
ax_zoom4.plot(x_exact, density_exact, label='Exact', linestyle=':', color='blue')
ax_zoom4.set_title("Shock Wave")
ax_zoom4.set_xlabel("x")
ax_zoom4.set_ylabel("Density")
ax_zoom4.set_xlim(0.846, 0.855)
ax_zoom4.set_ylim(0.1, 0.275)
ax_zoom4.legend()
ax_zoom4.spines['top'].set_visible(False)
ax_zoom4.spines['right'].set_visible(False)

# Adjust layout and save the figure
plt.tight_layout(rect=[0, 0, 1, 0.96])  # Leave space for the suptitle
plt.subplots_adjust(hspace=0.45, wspace=0.15)  # Reduce space between rows and columns
plt.savefig("Calibration-TENO")
plt.show()
