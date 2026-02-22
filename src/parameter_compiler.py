import re

# ==============================================================================
# 1. FIXED-POINT CONVERTER (The Bridge to Hardware)
# ==============================================================================
def to_q16_16(value):
    """Converts a float to a 32-bit Q16.16 hex string."""
    if value is None: return "00000000"
    
    # Handle constants (like periodicity n) that are integers
    if isinstance(value, int):
        return f"{value:x}" 
        
    # Fixed-Point Math: Multiply by 2^16 (65536)
    fixed_val = int(value * 65536)
    
    # Handle negative numbers (Two's Complement)
    if fixed_val < 0:
        fixed_val = (1 << 32) + fixed_val
        
    return f"{fixed_val:08X}" # Returns 8-char hex string (e.g., "00018000")

# ==============================================================================
# 2. CHARMM FILE PARSER
# ==============================================================================
class ForceField:
    def __init__(self):
        self.bonds = {}      # Key: (Type1, Type2) -> (kb, r0)
        self.angles = {}     # Key: (Type1, Type2, Type3) -> (ktheta, theta0)
        self.dihedrals = {}  # Key: (Type1, Type2, Type3, Type4) -> (kchi, n, delta)
        self.atom_types = {} # Key: AtomName -> (AtomType, Charge)

    def load_prm(self, filename):
        """Parses the .prm parameter file."""
        section = None
        with open(filename, 'r') as f:
            for line in f:
                # 1. Strip inline comments and whitespace safely
                line = line.split('!')[0].strip()
                if not line: continue
                
                # 2. Detect ALL possible CHARMM Sections to prevent bleed-over
                if line.startswith('BONDS'): section = 'BONDS'; continue
                if line.startswith('ANGLES'): section = 'ANGLES'; continue
                if line.startswith('DIHEDRALS'): section = 'DIHEDRALS'; continue
                if line.startswith('IMPROPER'): section = 'IMPROPER'; continue
                if line.startswith('CMAP'): section = 'CMAP'; continue
                if line.startswith('NONBONDED'): section = 'NONBONDED'; continue
                if line.startswith('NBFIX'): section = 'NBFIX'; continue
                if line.startswith('HBOND'): section = 'HBOND'; continue

                parts = line.split()
                
                # 3. Safe Parsing: Check column counts before extracting
                if section == 'BONDS' and len(parts) >= 4:
                    # Format: Atom1 Atom2 Kb r0
                    key = tuple(sorted((parts[0], parts[1]))) 
                    self.bonds[key] = (float(parts[2]), float(parts[3]))

                elif section == 'ANGLES' and len(parts) >= 5:
                    # Format: Atom1 Atom2 Atom3 Ktheta Theta0
                    key = (min(parts[0], parts[2]), parts[1], max(parts[0], parts[2]))
                    self.angles[key] = (float(parts[3]), float(parts[4]))

                elif section == 'DIHEDRALS' and len(parts) >= 7:
                    # Format: A B C D Kchi n delta
                    key = tuple(parts[0:4])
                    self.dihedrals[key] = (float(parts[4]), int(parts[5]), float(parts[6]))

    def load_rtf(self, filename):
        """Parses the .rtf topology file to get Atom Types and Charges."""
        current_residue = None
        with open(filename, 'r') as f:
            for line in f:
                parts = line.split()
                if not parts: continue
                
                if parts[0] == 'RESI':
                    current_residue = parts[1]
                
                elif parts[0] == 'ATOM':
                    # Format: ATOM Name Type Charge
                    # Map: (Residue, AtomName) -> (AtomType, Charge)
                    atom_name = parts[1]
                    atom_type = parts[2]
                    charge = float(parts[3])
                    self.atom_types[(current_residue, atom_name)] = (atom_type, charge)

    # ==========================================================================
    # 3. PARAMETER LOOKUP LOGIC
    # ==========================================================================
    def get_bond_params(self, t1, t2):
        key = tuple(sorted((t1, t2)))
        params = self.bonds.get(key)
        if params: return params # (kb, r0)
        return (0.0, 0.0) # Default if not found

    def get_angle_params(self, t1, t2, t3):
        # Center atom (t2) must match. Outer atoms (t1, t3) can swap.
        key = (min(t1, t3), t2, max(t1, t3))
        params = self.angles.get(key)
        if params: return params # (ktheta, theta0)
        return (0.0, 90.0) 

    def get_dihedral_params(self, t1, t2, t3, t4):
        # Try Exact Match
        key = (t1, t2, t3, t4)
        if key in self.dihedrals: return self.dihedrals[key]
        
        # Try Reverse Match
        key_rev = (t4, t3, t2, t1)
        if key_rev in self.dihedrals: return self.dihedrals[key_rev]
        
        # Try Wildcard (X, t2, t3, X) - Common in CHARMM
        key_wild = ('X', t2, t3, 'X')
        if key_wild in self.dihedrals: return self.dihedrals[key_wild]
        
        return (0.0, 1, 180.0) # Default: No barrier, Trans conformation

# ==============================================================================
# 4. COMPILER MAIN ROUTINE
# ==============================================================================
def compile_hex_file(ff, atom_list, output_filename="forcefield.hex"):
    """
    Generates the HEX file for parameter_ram.v
    atom_list: List of tuples [(ResName, AtomName), ...] representing the linear chain.
    """
    print(f"Compiling {len(atom_list)} atoms into {output_filename}...")
    
    with open(output_filename, 'w') as f:
        # Sliding Window of 4 Atoms (A, B, C, D)
        for i in range(len(atom_list) - 3):
            # 1. Identify the 4 atoms in the window
            a1 = atom_list[i]
            a2 = atom_list[i+1]
            a3 = atom_list[i+2]
            a4 = atom_list[i+3]
            
            # 2. Get their CHARMM Atom Types (e.g., 'CT1', 'NH1')
            t1, q1 = ff.atom_types[a1]
            t2, q2 = ff.atom_types[a2]
            t3, q3 = ff.atom_types[a3]
            t4, q4 = ff.atom_types[a4]
            
            # 3. Look up Physics Parameters
            kb, r0          = ff.get_bond_params(t1, t2)        # Bond A-B
            kth, th0        = ff.get_angle_params(t1, t2, t3)   # Angle A-B-C
            kphi, n, phi0   = ff.get_dihedral_params(t1,t2,t3,t4) # Dihedral A-B-C-D
            
            # 4. Convert to Verilog Q16.16 Hex
            # Note: CHARMM angles are in Degrees. Need Radians? 
            # Ideally yes, but let's assume your ASIC expects degrees for now 
            # OR multiply by (3.14159/180) here. Let's do Degrees -> Radians:
            th0_rad = th0 * (3.14159 / 180.0)
            phi0_rad = phi0 * (3.14159 / 180.0)

            hex_r0 = to_q16_16(r0)
            hex_kb = to_q16_16(kb)
            hex_th0 = to_q16_16(th0_rad)
            hex_kth = to_q16_16(kth)
            hex_phi0 = to_q16_16(phi0_rad)
            hex_kphi = to_q16_16(kphi)
            hex_n = f"{n:1x}" # 4-bit integer, not Q16.16
            hex_qa = to_q16_16(q1)
            hex_qd = to_q16_16(q4)
            
            # 5. Format the line (260 bits total width)
            # The order MUST match your Verilog concatenation from left (MSB) to right (LSB):
            # {r0, kb, theta0, k_theta, phi0, k_phi, n, q_a, q_d}
            
            # Note: n is a 4-bit integer, so we use :1x to format it as a single hex character
            hex_n = f"{n:1x}"
            
            # Combine into a single 65-character string (65 chars * 4 bits = 260 bits)
            hex_line = f"{hex_r0}{hex_kb}{hex_th0}{hex_kth}{hex_phi0}{hex_kphi}{hex_n}{hex_qa}{hex_qd}"
            
            f.write(f"// Atom {i}: {a1}-{a2}-{a3}-{a4} ({t1}-{t2}-{t3}-{t4})\n")
            f.write(f"{hex_line}\n")

# ==============================================================================
# 5. EXECUTION
# ==============================================================================
if __name__ == "__main__":
    # Initialize System
    ff = ForceField()
    
    # Load the files (Assuming you saved them locally)
    try:
        ff.load_rtf("top_all36_prot.rtf") 
        ff.load_prm("par_all36_prot.prm")
        print("Files loaded successfully.")
    except FileNotFoundError:
        print("Error: Ensure .rtf and .prm files are in the folder.")
        exit()

# --- DEFINE TEST MOLECULE (10-Atom Alanine Fragment) ---
# Ensure this sequence is long enough to generate 10 parameter rows
# Ensure this sequence is long enough to generate 10+ parameter rows
    test_sequence = [
        ('ALA', 'N'), ('ALA', 'HN'), ('ALA', 'CA'), ('ALA', 'HA'), 
        ('ALA', 'CB'), ('ALA', 'C'), ('ALA', 'O'), ('ALA', 'N'), 
        ('ALA', 'HN'), ('ALA', 'CA'), ('ALA', 'HA'), ('ALA', 'CB'), ('ALA', 'C')
    ]
    
    compile_hex_file(ff, test_sequence, "forcefield_init.hex")
    print("Done! 'forcefield_init.hex' generated.")