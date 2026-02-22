`default_nettype none

module atom_type_table (
    input  wire [15:0] type_id,     // From your Residue DB
    output reg  [31:0] sigma_q16,   // Atomic Radius (Angstroms) in Q16.16
    output reg  [31:0] epsilon_q16  // Energy Well Depth (kcal/mol) in Q16.16
);

    // Using standard CHARMM/AMBER force field values as referenced in the spec [cite: 29, 30]
    always @(*) begin
        case (type_id)
            // --- Carbon Types ---
            16'h0001: begin // C (Carbonyl)
                sigma_q16   = 32'h0001B333; // 1.700 A
                epsilon_q16 = 32'h00001C28; // 0.110 kcal/mol
            end
            16'h0009, 16'h000A: begin // CT1/CT2 (Aliphatic Carbon)
                sigma_q16   = 32'h0001E8F5; // 1.910 A
                epsilon_q16 = 32'h00000E14; // 0.055 kcal/mol
            end

            // --- Nitrogen Types ---
            16'h0010: begin // NH1 (Peptide Nitrogen)
                sigma_q16   = 32'h0001D1EB; // 1.820 A
                epsilon_q16 = 32'h00002B02; // 0.170 kcal/mol
            end

            // --- Oxygen Types ---
            16'h0014: begin // O (Carbonyl Oxygen)
                sigma_q16   = 32'h0001B333; // 1.700 A
                epsilon_q16 = 32'h00001EB8; // 0.120 kcal/mol
            end
            16'h0016: begin // OH1 (Hydroxyl Oxygen)
                sigma_q16   = 32'h00019EB8; // 1.620 A
                epsilon_q16 = 32'h000028F5; // 0.160 kcal/mol
            end

            // --- Sulfur Types ---
            16'h0017: begin // S (Sulfur)
                sigma_q16   = 32'h00020000; // 2.000 A
                epsilon_q16 = 32'h00003AE1; // 0.230 kcal/mol
            end

            default: begin
                sigma_q16   = 32'h00018000; // 1.500 A (Default)
                epsilon_q16 = 32'h00001999; // 0.100 (Default)
            end
        endcase
    end
endmodule