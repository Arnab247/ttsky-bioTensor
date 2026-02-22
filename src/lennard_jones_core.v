module lennard_jones_core (
    input clk,
    input rst_n,
    input signed [31:0] sigma_sq,        // Q16.16
    input signed [31:0] epsilon_x24,     // Q16.16
    input signed [31:0] r2_inv,          // Q16.16
    output reg signed [31:0] f_lj        // Q16.16
);

    // --- Math Pipeline Registers ---
    reg signed [31:0] sr2, sr4, sr6, sr12;
    reg signed [31:0] force_term, force_eps;

    // --- Delay Registers (Shift Registers) ---
    // These carry the original inputs forward so they align with the math stages
    reg signed [31:0] eps_d1, eps_d2, eps_d3, eps_d4, eps_d5;
    reg signed [31:0] r2_inv_d1, r2_inv_d2, r2_inv_d3, r2_inv_d4, r2_inv_d5, r2_inv_d6;
    reg signed [31:0] sr2_d2; // Needs to wait 1 cycle for Stage 3
    reg signed [31:0] sr6_d4; // Needs to wait 1 cycle for Stage 5

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            sr2 <= 0; sr4 <= 0; sr6 <= 0; sr12 <= 0;
            force_term <= 0; force_eps <= 0; f_lj <= 0;
            eps_d1 <= 0; eps_d2 <= 0; eps_d3 <= 0; eps_d4 <= 0; eps_d5 <= 0;
            r2_inv_d1 <= 0; r2_inv_d2 <= 0; r2_inv_d3 <= 0; r2_inv_d4 <= 0; r2_inv_d5 <= 0; r2_inv_d6 <= 0;
            sr2_d2 <= 0; sr6_d4 <= 0;
        end else begin
            // --- STAGE 1: Calculate (sigma/r)^2 ---
            sr2 <= (64'(sigma_sq) * r2_inv) >>> 16;
            eps_d1 <= epsilon_x24;
            r2_inv_d1 <= r2_inv;

            // --- STAGE 2: Calculate (sigma/r)^4 ---
            sr4 <= (64'(sr2) * sr2) >>> 16;
            sr2_d2 <= sr2;          // Save sr2 for the next stage
            eps_d2 <= eps_d1;
            r2_inv_d2 <= r2_inv_d1;

            // --- STAGE 3: Calculate (sigma/r)^6 ---
            sr6 <= (64'(sr4) * sr2_d2) >>> 16;
            eps_d3 <= eps_d2;
            r2_inv_d3 <= r2_inv_d2;

            // --- STAGE 4: Calculate (sigma/r)^12 ---
            sr12 <= (64'(sr6) * sr6) >>> 16;
            sr6_d4 <= sr6;          // Save sr6 for the subtraction stage
            eps_d4 <= eps_d3;
            r2_inv_d4 <= r2_inv_d3;

            // --- STAGE 5: The Subtraction ---
            force_term <= (sr12 << 1) - sr6_d4;
            eps_d5 <= eps_d4;
            r2_inv_d5 <= r2_inv_d4;

            // --- STAGE 6: Apply Epsilon ---
            force_eps <= (64'(force_term) * eps_d5) >>> 16;
            r2_inv_d6 <= r2_inv_d5;

            // --- STAGE 7: Apply Final 1/r^2 ---
            f_lj <= (64'(force_eps) * r2_inv_d6) >>> 16;
        end
    end
endmodule