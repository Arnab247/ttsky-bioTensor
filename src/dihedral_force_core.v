`default_nettype none

module dihedral_force_core (
    input  wire        clk,
    input  wire        rst_n,
    input  wire        start,
    
    // Coordinates
    input  wire signed [31:0] xa, ya, za, xb, yb, zb, xc, yc, zc, xd, yd, zd,
    
    // Torsion Parameters
    input  wire signed [31:0] phi0, k_phi,
    input  wire [3:0]         n_period,
    
    // Outputs
    output reg signed [31:0] phi_out,
    output reg signed [31:0] fax, fay, faz, fbx, fby, fbz, fcx, fcy, fcz, fdx, fdy, fdz,
    output reg         valid_out, busy
);

    // --- State Definitions ---
    localparam S_IDLE           = 0,
               S_VEC_GEN        = 1,
               S_CROSS_N1_PREP  = 2, S_CROSS_N1_LATCH = 3,
               S_CROSS_N2_PREP  = 4, S_CROSS_N2_LATCH = 5,
               S_NORM_B2_TRIG   = 6, S_NORM_B2_WAIT   = 7,
               S_NORM_N1_TRIG   = 8, S_NORM_N1_WAIT   = 9,
               S_NORM_N2_TRIG   = 10, S_NORM_N2_WAIT  = 11,
               S_CROSS_M1_PREP  = 12, S_CROSS_M1_LATCH = 13,
               S_DOT_XY         = 14,
               S_ATAN2_TRIG     = 15, S_ATAN2_WAIT    = 16,
               S_FORCE_CALC_1   = 17, // dPhi
               S_FORCE_CALC_2   = 18, // Torque and Scaling Factor Prep
               S_FORCE_CALC_3   = 19, // Calculate Fa/Fd Magnitude
               S_FORCE_CALC_4   = 20, // Final Projection
               S_OUT            = 21;

    reg [4:0] state;
    
    // Registers
    reg signed [31:0] b1x, b1y, b1z, b2x, b2y, b2z, b3x, b3y, b3z;
    reg signed [31:0] n1x, n1y, n1z, n2x, n2y, n2z, m1x, m1y, m1z, u_b2x, u_b2y, u_b2z;
    reg signed [31:0] dot_x, dot_y, d_phi, torque;
    
    // Scaling Registers
    reg signed [31:0] mag_sq_b2, mag_b2;
    reg signed [31:0] inv_mag_b2, inv_mag_n1, inv_mag_n2;
    reg signed [31:0] coeff_a, coeff_d;

    // --- Sub-modules ---
    
    // 1. Cross Product
    reg  signed [31:0] cp_v1x, cp_v1y, cp_v1z, cp_v2x, cp_v2y, cp_v2z;
    wire signed [31:0] cp_rx, cp_ry, cp_rz;
    vector_cross_product u_cp (.Ax(cp_v1x), .Ay(cp_v1y), .Az(cp_v1z), .Bx(cp_v2x), .By(cp_v2y), .Bz(cp_v2z), .Rx(cp_rx), .Ry(cp_ry), .Rz(cp_rz));

    // 2. Normalizer (Multiplexed)
    reg signed [31:0] norm_x, norm_y, norm_z;
    reg  n_start; wire n_valid, n_busy; 
    wire signed [31:0] n_nx, n_ny, n_nz, n_inv_mag; // Note: Ensure vector_normalizer has inv_mag output!
    
    vector_normalizer u_norm (
        .clk(clk), .rst_n(rst_n), .start(n_start), 
        .vx(norm_x), .vy(norm_y), .vz(norm_z),
        .nx(n_nx), .ny(n_ny), .nz(n_nz), 
        .inv_mag(n_inv_mag), // <--- CRITICAL CONNECTION
        .valid_out(n_valid), .busy(n_busy)
    );

    // 3. CORDIC Atan2
    reg  c_start; wire c_valid, c_busy; wire signed [31:0] c_theta;
    cordic_atan2 u_atan (.clk(clk), .rst_n(rst_n), .start(c_start), .x_in(dot_x), .y_in(dot_y), .theta_out(c_theta), .valid_out(c_valid), .busy(c_busy));

    function automatic signed [31:0] qmult(input signed [31:0] a, input signed [31:0] b);
        reg signed [63:0] temp;
        begin
            temp = $signed({{32{a[31]}}, a}) * $signed({{32{b[31]}}, b});
            qmult = temp[47:16];
        end
    endfunction

    // --- Main FSM ---
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            state <= S_IDLE; {valid_out, busy, n_start, c_start} <= 0;
        end else begin
            case (state)
                S_IDLE: begin
                    valid_out <= 0;
                    if (start) begin busy <= 1; state <= S_VEC_GEN; end
                end

                S_VEC_GEN: begin
                    b1x <= xb - xa; b1y <= yb - ya; b1z <= zb - za;
                    b2x <= xc - xb; b2y <= yc - yb; b2z <= zc - zb;
                    b3x <= xd - xc; b3y <= yd - yc; b3z <= zd - zc;
                    state <= S_CROSS_N1_PREP;
                end

                S_CROSS_N1_PREP: begin
                    cp_v1x <= b1x; cp_v1y <= b1y; cp_v1z <= b1z;
                    cp_v2x <= b2x; cp_v2y <= b2y; cp_v2z <= b2z;
                    state <= S_CROSS_N1_LATCH;
                end
                
                S_CROSS_N1_LATCH: begin
                    n1x <= cp_rx; n1y <= cp_ry; n1z <= cp_rz;
                    state <= S_CROSS_N2_PREP;
                end

                S_CROSS_N2_PREP: begin
                    cp_v1x <= b2x; cp_v1y <= b2y; cp_v1z <= b2z;
                    cp_v2x <= b3x; cp_v2y <= b3y; cp_v2z <= b3z;
                    state <= S_CROSS_N2_LATCH;
                end

                S_CROSS_N2_LATCH: begin
                    n2x <= cp_rx; n2y <= cp_ry; n2z <= cp_rz;
                    // Prepare Multiplexing: Send B2 first
                    norm_x <= b2x; norm_y <= b2y; norm_z <= b2z;
                    state <= S_NORM_B2_TRIG;
                end

                // --- Normalization 1: B2 ---
                S_NORM_B2_TRIG: begin n_start <= 1; state <= S_NORM_B2_WAIT; end
                S_NORM_B2_WAIT: begin
                    n_start <= 0;
                    if (n_valid) begin
                        u_b2x <= n_nx; u_b2y <= n_ny; u_b2z <= n_nz;
                        inv_mag_b2 <= n_inv_mag; // Capture 1/|B2|
                        
                        // Calculate |B2| magnitude for scaling: |B2| = (B2 dot B2) * (1/|B2|)
                        // Note: doing dot product here inline
                        mag_sq_b2 <= qmult(b2x, b2x) + qmult(b2y, b2y) + qmult(b2z, b2z);
                        
                        norm_x <= n1x; norm_y <= n1y; norm_z <= n1z;
                        state <= S_NORM_N1_TRIG;
                    end
                end

                // --- Normalization 2: N1 ---
                S_NORM_N1_TRIG: begin n_start <= 1; state <= S_NORM_N1_WAIT; end
                S_NORM_N1_WAIT: begin
                    n_start <= 0;
                    if (n_valid) begin
                         if (n1x == 0 && n1y == 0 && n1z == 0) begin
                            n1x <= 32'h00010000; n1y <= 0; n1z <= 0;
                        end else begin
                            n1x <= n_nx; n1y <= n_ny; n1z <= n_nz; // N1 is now unit vector
                        end
                        inv_mag_n1 <= n_inv_mag; // Capture 1/|N1|
                        
                        norm_x <= n2x; norm_y <= n2y; norm_z <= n2z;
                        state <= S_NORM_N2_TRIG;
                    end
                end

                // --- Normalization 3: N2 ---
                S_NORM_N2_TRIG: begin n_start <= 1; state <= S_NORM_N2_WAIT; end
                S_NORM_N2_WAIT: begin
                    n_start <= 0;
                    if (n_valid) begin
                        if (n2x == 0 && n2y == 0 && n2z == 0) begin
                            n2x <= 32'h00010000; n2y <= 0; n2z <= 0;
                        end else begin
                            n2x <= n_nx; n2y <= n_ny; n2z <= n_nz; // N2 is now unit vector
                        end
                        inv_mag_n2 <= n_inv_mag; // Capture 1/|N2|
                        state <= S_CROSS_M1_PREP;
                    end
                end

                S_CROSS_M1_PREP: begin
                    cp_v1x <= n1x; cp_v1y <= n1y; cp_v1z <= n1z;
                    cp_v2x <= u_b2x; cp_v2y <= u_b2y; cp_v2z <= u_b2z;
                    state <= S_CROSS_M1_LATCH;
                end

                S_CROSS_M1_LATCH: begin
                    m1x <= cp_rx; m1y <= cp_ry; m1z <= cp_rz;
                    state <= S_DOT_XY;
                end

                S_DOT_XY: begin
                    dot_x <= qmult(n1x, n2x) + qmult(n1y, n2y) + qmult(n1z, n2z);
                    dot_y <= qmult(m1x, n2x) + qmult(m1y, n2y) + qmult(m1z, n2z);
                    state <= S_ATAN2_TRIG;
                end

                S_ATAN2_TRIG: begin c_start <= 1; state <= S_ATAN2_WAIT; end
                S_ATAN2_WAIT: begin
                    c_start <= 0;
                    if (c_valid) begin
                        phi_out <= c_theta;
                        state   <= S_FORCE_CALC_1;
                    end
                end

                // --- PHYSICALLY ACCURATE FORCE PIPELINE ---
                S_FORCE_CALC_1: begin
                    d_phi <= c_theta - phi0;
                    // Calculate Real |B2| = |B2|^2 * (1/|B2|)
                    mag_b2 <= qmult(mag_sq_b2, inv_mag_b2);
                    state <= S_FORCE_CALC_2;
                end

                S_FORCE_CALC_2: begin
                    torque <= qmult(k_phi, d_phi);
                    state <= S_FORCE_CALC_3;
                end

                S_FORCE_CALC_3: begin
                    // Wilson-Decius Scaling:
                    // Force = Torque * (|Bond_BC| / |Normal|)
                    // We use inv_mag_n1 which is (1 / |N1_unnormalized|)
                    
                    // Coefficient A = Torque * |B2| * (1/|N1|)
                    coeff_a <= qmult(torque, qmult(mag_b2, inv_mag_n1));
                    
                    // Coefficient D = -Torque * |B2| * (1/|N2|)
                    coeff_d <= -qmult(torque, qmult(mag_b2, inv_mag_n2));
                    
                    state <= S_FORCE_CALC_4;
                end

                S_FORCE_CALC_4: begin
                    // Apply Coefficients to Unit Normals
                    fax <= qmult(coeff_a, n1x); 
                    fay <= qmult(coeff_a, n1y); 
                    faz <= qmult(coeff_a, n1z);

                    fdx <= qmult(coeff_d, n2x); 
                    fdy <= qmult(coeff_d, n2y); 
                    fdz <= qmult(coeff_d, n2z);

                    // Simplified Conservation (Forces on B/C oppose A/D)
                    // Note: Full angular momentum conservation requires 
                    // complex coefficients for Fb/Fc involving cotangents.
                    // This linear approx prevents drift but may leave residual torque.
                    fbx <= -qmult(coeff_a, n1x); 
                    fby <= -qmult(coeff_a, n1y); 
                    fbz <= -qmult(coeff_a, n1z);
                    
                    fcx <= -qmult(coeff_d, n2x); 
                    fcy <= -qmult(coeff_d, n2y); 
                    fcz <= -qmult(coeff_d, n2z);

                    state <= S_OUT;
                end

                S_OUT: begin
                    valid_out <= 1; busy <= 0; state <= S_IDLE;
                end
            endcase
        end
    end
endmodule