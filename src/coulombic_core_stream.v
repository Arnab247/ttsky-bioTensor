`default_nettype none

module coulombic_core_stream #(
    parameter signed [31:0] KC = 32'h014C1000 // 332.06 in Q16.16
)(
    input  wire clk,
    input  wire rst_n,
    input  wire signed [31:0] q_i,      // Q16.16
    input  wire signed [31:0] q_j,      // Q16.16
    input  wire signed [31:0] r2_inv,   // Q16.16 (From shared NBPipe)
    output reg  signed [31:0] f_scalar  // Q16.16
);

    // --- Pipeline Registers ---
    reg signed [31:0] stage1_qq;
    reg signed [31:0] stage1_r2_inv; // NEW: Delay register for alignment
    reg signed [31:0] stage2_raw_force;

    // --- Verilog-2001 Safe Q16.16 Multiplier ---
    function automatic signed [31:0] qmult(input signed [31:0] a, input signed [31:0] b);
        reg signed [63:0] temp;
        begin
            temp = $signed({{32{a[31]}}, a}) * $signed({{32{b[31]}}, b});
            qmult = temp[47:16];
        end
    endfunction

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            stage1_qq        <= 0;
            stage1_r2_inv    <= 0;
            stage2_raw_force <= 0;
            f_scalar         <= 0;
        end else begin
            // STAGE 1: Charge Product & Delay r2_inv
            stage1_qq     <= qmult(q_i, q_j);
            stage1_r2_inv <= r2_inv; // Hold r2_inv for 1 cycle

            // STAGE 2: Geometric Decay 
            // Now stage1_qq and stage1_r2_inv are perfectly aligned!
            stage2_raw_force <= qmult(stage1_qq, stage1_r2_inv);

            // STAGE 3: Unit Scaling
            f_scalar <= qmult(stage2_raw_force, KC);
        end
    end

endmodule