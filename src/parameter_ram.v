module parameter_ram (
    input wire clk,
    input wire [9:0] addr, // Matches your scan_idx width (supports up to 1024 atoms)

    // 2-Body Bond Parameters
    output reg signed [31:0] r0_out,       // Target Bond Length
    output reg signed [31:0] kb_out,       // Bond Stiffness

    // 3-Body Angle Parameters
    output reg signed [31:0] theta0_out,   // Target Angle
    output reg signed [31:0] k_theta_out,  // Angle Stiffness

    // 4-Body Dihedral Parameters
    output reg signed [31:0] phi0_out,     // Target Torsion Phase
    output reg signed [31:0] k_phi_out,    // Torsion Stiffness
    output reg [3:0]         n_period_out, // Periodicity (1, 2, 3, etc.)

    // Electrostatics (Charges)
    output reg signed [31:0] q_a_out,      // Partial Charge of Atom A
    output reg signed [31:0] q_d_out       // Partial Charge of Atom D
);

    // ---------------------------------------------------------
    // WIDE MEMORY ARRAY
    // ---------------------------------------------------------
    // Total Width: 32*8 + 4 = 260 bits wide
    // Depth: 1024 entries (one set of params per atom window)
    reg [259:0] mem [0:1023]; 

    // ---------------------------------------------------------
    // SYNCHRONOUS READ
    // ---------------------------------------------------------
    always @(posedge clk) begin
        // Fetch the massive 260-bit word and slice it into wires
        {r0_out, kb_out, theta0_out, k_theta_out, phi0_out, k_phi_out, n_period_out, q_a_out, q_d_out} <= mem[addr];
    end

    // ---------------------------------------------------------
    // SIMULATION INITIALIZATION (The "Stub" Data)
    // ---------------------------------------------------------
    // This pre-loads the memory with standard Carbon-Carbon backbone
    // values so your simulation works immediately without a host.
    integer i;
    initial begin
        // 1. Initialize EVERYTHING to safe zero-energy defaults first
        for (i = 0; i < 1024; i = i + 1) begin
            mem[i] = 260'd0; 
        end
        
        // 2. Load the compiled hex file
        $readmemh("forcefield_init.hex", mem);
    end

endmodule