`default_nettype none

module non_bonded_pipeline (
    input  wire        clk,
    input  wire        rst_n,
    input  wire        valid_in,   
    
    input  wire signed [31:0] xi, yi, zi,
    input  wire signed [31:0] xj, yj, zj,
    
    input  wire signed [31:0] q_i, q_j,
    input  wire signed [31:0] sigma_sq, eps_x24,
    
    output wire                valid_out, 
    output wire  signed [31:0] fx_nb, fy_nb, fz_nb 
);

    // STAGE 0: Input capture (ensure coordinates are stable before calculation)
    reg signed [31:0] xi_r, yi_r, zi_r;
    reg signed [31:0] xj_r, yj_r, zj_r;
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            xi_r <= 0; yi_r <= 0; zi_r <= 0;
            xj_r <= 0; yj_r <= 0; zj_r <= 0;
        end else begin
            xi_r <= xi; yi_r <= yi; zi_r <= zi;
            xj_r <= xj; yj_r <= yj; zj_r <= zj;
        end
    end
    
    // STAGE 1: Distance Calculation
    reg signed [31:0] dx, dy, dz;
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin dx <= 0; dy <= 0; dz <= 0; end
        else begin dx <= xi_r - xj_r; dy <= yi_r - yj_r; dz <= zi_r - zj_r; end
    end

    // STAGE 2: R2 
    wire signed [63:0] dx2 = $signed(dx) * $signed(dx);
    wire signed [63:0] dy2 = $signed(dy) * $signed(dy);
    wire signed [63:0] dz2 = $signed(dz) * $signed(dz);
    
    reg signed [31:0] r2;
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) r2 <= 0;
        else r2 <= $signed((dx2 + dy2 + dz2) >>> 16);  // Always update!
    end

    // ==============================================================
    // STAGE 1-8: Parameters (Stretched for 6-Cycle Wrapper)
    // ==============================================================
    reg signed [31:0] qi_sr [0:8], qj_sr [0:8], sig_sr [0:8], eps_sr [0:8];
    integer j;
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) for (j=0; j<=8; j=j+1) begin qi_sr[j] <= 0; qj_sr[j] <= 0; sig_sr[j] <= 0; eps_sr[j] <= 0; end
        else begin
            qi_sr[0] <= q_i; qj_sr[0] <= q_j; sig_sr[0] <= sigma_sq; eps_sr[0] <= eps_x24;
            for (j = 1; j <= 8; j = j + 1) begin
                qi_sr[j] <= qi_sr[j-1]; qj_sr[j] <= qj_sr[j-1];
                sig_sr[j] <= sig_sr[j-1]; eps_sr[j] <= eps_sr[j-1];
            end
        end
    end

    // STAGE 3-8: Reciprocal Wrapper (6 cycles)
    wire [31:0] r2_inv;
    wire recip_valid_out;  
    reciprocal_wrapper u_shared_inv (
        .clk(clk), 
        .rst_n(rst_n), 
        .valid_in(v_sr[2]),  // Starts at Stage 2
        .x_in(r2), 
        .y_out(r2_inv),      // Finishes at Stage 8
        .valid_out(recip_valid_out)
    );

    // --- PARALLEL MATH ENGINES (Start at Stage 2) ---
    // Your existing reciprocal unit:
    // reciprocal_wrapper u_recip ( ... );

    // NEW: The 1/r engine
    wire signed [31:0] inv_r_out;
    
    // The new 1/r engine (Starts at Stage 2, Finishes at Stage 8)
    inv_sqrt_direct u_inv_sqrt (
        .clk(clk),
        .rst_n(rst_n),
        .valid_in(v_sr[2]), 
        .x(r2),             
        .y_out(inv_r_out),
        .valid_out()        
    );


    // Shift register to carry 1/r from Stage 8 down to the end of the pipeline
    reg signed [31:0] inv_r_sr [8:16]; 
    // Shift 1/r down the pipeline alongside the vectors
    integer ir_idx;
    always @(posedge clk) begin
        inv_r_sr[8] <= inv_r_out;
        for (ir_idx = 9; ir_idx <= 16; ir_idx = ir_idx + 1) begin
            inv_r_sr[ir_idx] <= inv_r_sr[ir_idx - 1];
        end
    end


    // ==============================================================
    // STAGE 8-15: Core Physics (Using Stage 8 Parameters)
    // ==============================================================
    wire signed [31:0] f_coulomb_raw, f_lj;
    coulombic_core_stream u_coulomb (.clk(clk), .rst_n(rst_n), .q_i(qi_sr[8]), .q_j(qj_sr[8]), .r2_inv(r2_inv), .f_scalar(f_coulomb_raw));
    lennard_jones_core u_lj (.clk(clk), .rst_n(rst_n), .sigma_sq(sig_sr[8]), .epsilon_x24(eps_sr[8]), .r2_inv(r2_inv), .f_lj(f_lj));

    // STAGE 15: Final Scalar Summation
    reg signed [31:0] f_scalar_reg;
    reg signed [31:0] f_c_d1, f_c_d2, f_c_d3, f_c_d4;
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            f_c_d1 <= 0; f_c_d2 <= 0; f_c_d3 <= 0; f_c_d4 <= 0;
            f_scalar_reg <= 0;
        end else begin
            f_c_d1 <= f_coulomb_raw; f_c_d2 <= f_c_d1; f_c_d3 <= f_c_d2; f_c_d4 <= f_c_d3;
            f_scalar_reg <= f_c_d4 + f_lj; // Valid at Stage 16
        end
    end

 // ==============================================================
    // STAGE 1-18: Vector Shadowing (Stretched to match new pipeline)
    // ==============================================================
    reg signed [31:0] dx_sr [0:17], dy_sr [0:17], dz_sr [0:17];
    integer k;
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            for (k=0; k<=17; k=k+1) begin dx_sr[k] <= 0; dy_sr[k] <= 0; dz_sr[k] <= 0; end
        end else begin
            dx_sr[0] <= dx; dy_sr[0] <= dy; dz_sr[0] <= dz;
            for (k = 1; k <= 17; k = k + 1) begin // FIXED: Now shifts all the way to 17
                dx_sr[k] <= dx_sr[k-1]; dy_sr[k] <= dy_sr[k-1]; dz_sr[k] <= dz_sr[k-1];
            end
        end
    end

    // --- 18-CYCLE CONTROL PATH ---
    reg [17:0] v_sr; 
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) v_sr <= 18'b0; 
        else v_sr <= {v_sr[16:0], valid_in};
    end
    
    // ==============================================================
    // STAGE 17-18: Final Vector Force Calculation 
    // ==============================================================
    function automatic signed [31:0] qmult(input signed [31:0] a, input signed [31:0] b);
        reg signed [63:0] temp;
        begin
            temp = $signed({{32{a[31]}}, a}) * $signed({{32{b[31]}}, b});
            qmult = temp[47:16];
        end
    endfunction

    // NEW STAGE 17: Normalize the force (F_scalar * 1/r)
    reg signed [31:0] f_norm_reg;
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            f_norm_reg <= 0;
        end else begin
            // f_scalar_reg perfectly aligns with inv_r_sr[16] at Stage 17
            f_norm_reg <= qmult(f_scalar_reg, inv_r_sr[16]); 
        end
    end
    
    // FINAL OUTPUT (Stage 18): Project normalized force into 3D space
    // f_norm_reg is ready at Stage 18. dx_sr[15] perfectly aligns with Stage 18.
    assign fx_nb = -qmult(f_norm_reg, dx_sr[15]);
    assign fy_nb = -qmult(f_norm_reg, dy_sr[15]);
    assign fz_nb = -qmult(f_norm_reg, dz_sr[15]);
    
    // Output valid signal perfectly aligned with the math (Stage 18)
    assign valid_out = v_sr[17];

endmodule