import matplotlib.pyplot as plt
import re

def parse_coordinates(filename):
    initial_coords = []
    final_coords = []
    current_list = None

    with open(filename, 'r') as f:
        for line in f:
            if "INITIAL COORDINATES" in line:
                current_list = initial_coords
            elif "FINAL COORDINATES" in line:
                current_list = final_coords
            elif "Atom" in line and current_list is not None:
                # Regex to extract X, Y, Z floating point numbers
                match = re.search(r'X=\s*([-\d.]+),\s*Y=\s*([-\d.]+),\s*Z=\s*([-\d.]+)', line)
                if match:
                    x, y, z = map(float, match.groups())
                    current_list.append((x, y, z))
                    
    return initial_coords, final_coords

# Read the data from your Verilog log
initial, final = parse_coordinates('sim_output2.txt')

if not initial or not final:
    print("Could not find coordinate data. Make sure you saved the terminal output to 'sim_output.txt'!")
    exit()

# Setup 3D Plot
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

# Unpack coordinates
ix, iy, iz = zip(*initial)
fx, fy, fz = zip(*final)

# Plot Initial (Red Dots) and Final (Blue Stars)
ax.scatter(ix, iy, iz, c='red', s=50, label='Initial (T=0)', alpha=0.5)
ax.scatter(fx, fy, fz, c='blue', marker='o', s=80, label='Final (T=5)')

# ==========================================
# NEW: DRAW THE CHEMICAL BONDS
# ==========================================
# List of (Atom A, Atom B) indices that share a chemical bond
bonds = [
    (0, 1), # N  - HN
    (0, 2), # N  - CA
    (2, 3), # CA - HA
    (2, 4), # CA - CB
    (2, 5), # CA - C
    (5, 6), # C  - O
    (5, 7), # C  - N  (The Peptide Bond bridging the two amino acids!)
    (7, 8), # N  - HN
    (7, 9)  # N  - CA
]

# Draw Initial Bonds (Faint Red)
for (atom1, atom2) in bonds:
    ax.plot([ix[atom1], ix[atom2]], 
            [iy[atom1], iy[atom2]], 
            [iz[atom1], iz[atom2]], 
            color='red', alpha=0.3, linewidth=2)

# Draw Final Bonds (Solid Blue)
for (atom1, atom2) in bonds:
    ax.plot([fx[atom1], fx[atom2]], 
            [fy[atom1], fy[atom2]], 
            [fz[atom1], fz[atom2]], 
            color='blue', alpha=0.8, linewidth=3)
# ==========================================

# Draw arrows showing the displacement path of each atom (optional, faint)
for i in range(len(initial)):
    ax.plot([ix[i], fx[i]], [iy[i], fy[i]], [iz[i], fz[i]], 'k:', alpha=0.3)
    # Add atom labels slightly offset from the final positions
    ax.text(fx[i], fy[i], fz[i], f"  A{i}", size=9, zorder=1, color='k')

ax.set_title("Hardware Minimization: Alanine Dipeptide Structure")
ax.set_xlabel("X (Å)")
ax.set_ylabel("Y (Å)")
ax.set_zlabel("Z (Å)")
ax.legend()
plt.show()