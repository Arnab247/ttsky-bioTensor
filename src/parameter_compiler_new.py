import re
import math

# ==============================================================================
# 1. FIXED-POINT CONVERTER 
# ==============================================================================
def to_q16_16(value):
    """Converts a float to a 32-bit Q16.16 hex string."""
    if value is None: return "00000000"
    if isinstance(value, int): return f"{value:x}" 
        
    fixed_val = int(value * 65536)
    if fixed_val < 0:
        fixed_val = (1 << 32) + fixed_val
        
    return f"{fixed_val:08X}"

# ==============================================================================
# 2. CHARMM FILE PARSER
# ==============================================================================
class ForceField:
    def __init__(self):
        self.bonds = {}      
        self.angles = {}     
        self.dihedrals = {}  
        self.atom_types = {} 
        self.nonbonded = {}  # NEW: Key: AtomType -> (Epsilon, Rmin/2)

    def load_prm(self, filename):
        """Parses the .prm parameter file."""
        section = None
        with open(filename, 'r') as f:
            for line in f:
                line = line.split('!')[0].strip()
                if not line: continue
                
                # Detect Sections
                if line.startswith('BONDS'): section = 'BONDS'; continue
                if line.startswith('ANGLES'): section = 'ANGLES'; continue
                if line.startswith('DIHEDRALS'): section = 'DIHEDRALS'; continue
                if line.startswith('NONBONDED'): section = 'NONBONDED'; continue

                parts = line.split()
                
                # --- NEW: Parse Non-Bonded Parameters ---
                if section == 'NONBONDED' and len(parts) >= 4:
                    # CHARMM format: AtomType  Ignored  Epsilon  Rmin/2
                    # The try-except block ignores header lines containing words like 'ctofnb'
                    try:
                        atom_type = parts[0]
                        epsilon = abs(float(parts[2])) 
                        rmin_half = float(parts[3])
                        self.nonbonded[atom_type] = (epsilon, rmin_half)
                    except ValueError:
                        # If parts[2] isn't a number, skip this line and keep going
                        continue

    def load_rtf(self, filename):
        """Parses the .rtf topology file to get Atom Types and Charges."""
        current_residue = None
        with open(filename, 'r') as f:
            for line in f:
                line = line.split('!')[0].strip() # Clean comments
                parts = line.split()
                if not parts: continue
                
                if parts[0] == 'RESI':
                    current_residue = parts[1]
                
                elif parts[0] == 'ATOM':
                    atom_name = parts[1]
                    atom_type = parts[2]
                    charge = float(parts[3])
                    self.atom_types[(current_residue, atom_name)] = (atom_type, charge)

    # ==========================================================================
    # 3. NON-BONDED LOOKUP & MATH CONVERSION
    # ==========================================================================
    def get_nonbonded_hardware_params(self, res_name, atom_name):
        """Fetches and converts parameters to hardware-ready formats."""
        # 1. Get Type and Charge
        atom_type, q = self.atom_types.get((res_name, atom_name), ("UNKNOWN", 0.0))
        
        # 2. Get VDW params
        epsilon, rmin_half = self.nonbonded.get(atom_type, (0.0, 0.0))
        
        # 3. Math Conversions for Hardware
        # Sigma = Rmin / (2^(1/6)) = (2 * Rmin_half) / 1.122462
        sigma = (2.0 * rmin_half) / 1.12246204831
        sigma_sq = sigma * sigma
        
        # Epsilon scaled by 24 for the derivative calculation in hardware
        eps_x24 = epsilon * 24.0
        
        return q, sigma_sq, eps_x24, atom_type

# ==============================================================================
# 4. COMPILER MAIN ROUTINE (Non-Bonded LUT Generation)
# ==============================================================================
def compile_nonbonded_lut(ff, atom_list, output_filename="nonbonded_lut.hex"):
    """
    Generates a memory file for the Parameter LUT.
    Each line corresponds to an atom's static parameters: {Q, Sigma^2, 24*Epsilon}
    """
    print(f"Compiling Non-Bonded parameters for {len(atom_list)} atoms into {output_filename}...")
    
    with open(output_filename, 'w') as f:
        for i, (res_name, atom_name) in enumerate(atom_list):
            
            # Fetch and convert math
            q, sigma_sq, eps_x24, atom_type = ff.get_nonbonded_hardware_params(res_name, atom_name)
            
            # Convert to Q16.16 Hex Strings
            hex_q = to_q16_16(q)
            hex_sig_sq = to_q16_16(sigma_sq)
            hex_eps24 = to_q16_16(eps_x24)
            
            # Write to file (96 bits total: 32 bit Q, 32 bit Sig^2, 32 bit Eps24)
            f.write(f"// Atom {i}: {res_name}-{atom_name} (Type: {atom_type}) | Q={q}, Sig^2={sigma_sq:.3f}, Eps*24={eps_x24:.3f}\n")
            f.write(f"{hex_q}{hex_sig_sq}{hex_eps24}\n")

# ==============================================================================
# 5. EXECUTION
# ==============================================================================
if __name__ == "__main__":
    ff = ForceField()
    
    try:
        ff.load_rtf("top_all36_prot.rtf") 
        ff.load_prm("par_all36_prot.prm")
        print("Force Field Files loaded successfully.")
    except FileNotFoundError:
        print("Error: Ensure .rtf and .prm files are in the folder.")
        exit()

    # --- DEFINE TEST MOLECULE (Alanine Dipeptide subset) ---
    test_sequence = [
        ('ALA', 'N'),
        ('ALA', 'CA'),
        ('ALA', 'CB'),  # Added sidechain carbon
        ('ALA', 'C'),
        ('ALA', 'O')
    ]
    
    compile_nonbonded_lut(ff, test_sequence, "nonbonded_lut.hex")
    print("Done! 'nonbonded_lut.hex' generated.")