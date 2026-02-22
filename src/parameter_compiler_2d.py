import math

# ==============================================================================
# 1. FIXED-POINT CONVERTER (Q16.16)
# ==============================================================================
def to_q16_16(value):
    """Converts a float to a 32-bit Q16.16 hex string."""
    if value is None: return "00000000"
    fixed_val = int(value * 65536)
    if fixed_val < 0:
        fixed_val = (1 << 32) + fixed_val
    return f"{fixed_val:08X}"

# ==============================================================================
# 2. FORCE FIELD PARSER
# ==============================================================================
class ForceField:
    def __init__(self):
        self.atom_types = {} 
        self.nonbonded = {}  

    def load_rtf(self, filename):
        current_res = None
        try:
            with open(filename, 'r') as f:
                for line in f:
                    line = line.split('!')[0].strip()
                    parts = line.split()
                    if not parts: continue
                    if parts[0] == 'RESI': current_res = parts[1]
                    elif parts[0] == 'ATOM' and current_res:
                        self.atom_types[(current_res, parts[1])] = (parts[2], float(parts[3]))
        except FileNotFoundError:
            print(f"Warning: {filename} not found.")

    def load_prm(self, filename):
        section = None
        try:
            with open(filename, 'r') as f:
                for line in f:
                    line = line.split('!')[0].strip()
                    if line.startswith('NONBONDED'): section = 'NONBONDED'; continue
                    parts = line.split()
                    if section == 'NONBONDED' and len(parts) >= 4:
                        try:
                            self.nonbonded[parts[0]] = (abs(float(parts[2])), float(parts[3]))
                        except ValueError: continue
        except FileNotFoundError:
            print(f"Warning: {filename} not found.")

# ==============================================================================
# 3. ASSET GENERATION (1D Identity & 2D Mixing Matrix)
# ==============================================================================
def compile_hardware_assets(ff, atom_list):
    # 1. Identify Unique Atom Types
    unique_types = []
    for res, name in atom_list:
        t_str, _ = ff.atom_types.get((res, name), ("UNKNOWN", 0.0))
        if t_str not in unique_types:
            unique_types.append(t_str)
    
    type_to_id = {t: i for i, t in enumerate(unique_types)}
    num_types = len(unique_types)
    
    # 2. Generate atom_identity.hex (1D)
    print(f"Generating identity hex for {len(atom_list)} atoms...")
    with open("atom_identity.hex", "w") as f:
        for res, name in atom_list:
            t_str, q = ff.atom_types.get((res, name), ("UNKNOWN", 0.0))
            tid = type_to_id[t_str]
            f.write(f"{to_q16_16(q)}{tid:02X}\n")

    # 3. Generate mixing_matrix.hex (2D)
    print(f"Generating 2D mixing matrix for {num_types} unique types...")
    with open("mixing_matrix.hex", "w") as f:
        for i in range(num_types):
            for j in range(num_types):
                eps_i, rmin_i = ff.nonbonded.get(unique_types[i], (0.0, 0.0))
                eps_j, rmin_j = ff.nonbonded.get(unique_types[j], (0.0, 0.0))
                
                # Lorentz-Berthelot Rules
                mixed_eps = math.sqrt(eps_i * eps_j)
                mixed_rmin = rmin_i + rmin_j
                mixed_sigma = mixed_rmin * 0.8908987 # Rmin * 2^(-1/6)
                
                sig_sq = mixed_sigma**2
                eps_24 = mixed_eps * 24.0
                
                f.write(f"{to_q16_16(sig_sq)}{to_q16_16(eps_24)}\n")

# ==============================================================================
# 4. EXECUTION (10-Atom Sequence)
# ==============================================================================
if __name__ == "__main__":
    ff = ForceField()
    # Replace with path to standard CHARMM force field files
    ff.load_rtf("top_all36_prot.rtf")
    ff.load_prm("par_all36_prot.prm")

    # 10-atom sequence representing two Alanine residues
    test_sequence = [
        ('ALA', 'N'),  ('ALA', 'HN'), ('ALA', 'CA'), ('ALA', 'HA'), ('ALA', 'CB'),
        ('ALA', 'C'),  ('ALA', 'O'),  ('ALA', 'N'),  ('ALA', 'HN'), ('ALA', 'CA')
    ]
    
    compile_hardware_assets(ff, test_sequence)
    print("Done. Generated atom_identity.hex and mixing_matrix.hex.")