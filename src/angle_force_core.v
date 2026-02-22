`default_nettype none

module angle_force_core (
    input  wire        clk,
    input  wire        rst_n,
    input  wire        start,
    input  wire signed [31:0] xa, ya, za, xb, yb, zb, xc, yc, zc,
    input  wire signed [31:0] theta0, k_theta,
    output reg signed [31:0] fax, fay, faz, fcx, fcy, fcz,
    output reg signed [31:0] fbx, fby, fbz, // Added Vertex Force
    output reg         valid_out,
    output reg         busy
);

    // --- State Definitions ---
    localparam S_IDLE = 0, S_VEC_GEN = 1;
    localparam S_NORM_BA = 2, S_NORM_BA_W = 3;
    localparam S_NORM_BC = 4, S_NORM_BC_W = 5;
    localparam S_ACOS = 6, S_ACOS_W = 7;
    localparam S_CROSS_1 = 8, S_CROSS_A = 9, S_CROSS_C = 10, S_OUT = 11;

    reg [3:0] state;
    reg signed [31:0] d_theta, uBAx, uBAy, uBAz, uBCx, uBCy, uBCz, nx, ny, nz;
    
    // Sub-module Handshakes
    reg signed [31:0] n_vx, n_vy, n_vz;
    wire signed [31:0] n_nx, n_ny, n_nz;
    reg n_start;
    wire n_valid, n_busy;

    reg signed [31:0] a_in;
    wire signed [31:0] a_out;
    reg a_start;
    wire a_valid, a_busy;

    reg signed [31:0] cp_v1x, cp_v1y, cp_v1z, cp_v2x, cp_v2y, cp_v2z;
    wire signed [31:0] cp_rx, cp_ry, cp_rz;

    reg signed [31:0] inv_len_ba, inv_len_bc;
    wire signed [31:0] n_inv_len;

    // Sub-module Instances
    acos_poly u_acos (.clk(clk), .rst_n(rst_n), .start(a_start), .x_in(a_in), .theta_out(a_out), .valid_out(a_valid), .busy(a_busy));
    vector_normalizer u_norm (.clk(clk), .rst_n(rst_n), .start(n_start), .vx(n_vx), .vy(n_vy), .vz(n_vz), .nx(n_nx), .ny(n_ny), .nz(n_nz), .inv_mag(n_inv_len), .valid_out(n_valid), .busy(n_busy));
    vector_cross_product u_cp (.Ax(cp_v1x), .Ay(cp_v1y), .Az(cp_v1z), .Bx(cp_v2x), .By(cp_v2y), .Bz(cp_v2z), .Rx(cp_rx), .Ry(cp_ry), .Rz(cp_rz));

    function automatic signed [31:0] qmult(input signed [31:0] a, input signed [31:0] b);
        reg signed [63:0] temp;
        begin
            temp = $signed({{32{a[31]}}, a}) * $signed({{32{b[31]}}, b});
            qmult = temp[47:16];
        end
    endfunction

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            state <= S_IDLE; {n_start, a_start, valid_out, busy} <= 0;
        end else begin
            case (state)
                S_IDLE: begin
                    valid_out <= 0; // Ensure pulse is only 1 cycle long
                    if (start) begin 
                        busy <= 1; 
                        state <= S_VEC_GEN; 
                    end
                end

                S_VEC_GEN: begin
                    n_vx <= xa - xb; n_vy <= ya - yb; n_vz <= za - zb;
                    state <= S_NORM_BA;
                end

                S_NORM_BA: begin
                    n_start <= 1;
                    if (n_busy) begin n_start <= 0; state <= S_NORM_BA_W; end
                end

                S_NORM_BA_W: if (n_valid) begin
                    {uBAx, uBAy, uBAz} <= {n_nx, n_ny, n_nz};
                    inv_len_ba <= n_inv_len; // Latch 1/|BA|
                    n_vx <= xc - xb; n_vy <= yc - yb; n_vz <= zc - zb;
                    state <= S_NORM_BC;
                end

                S_NORM_BC: begin
                    n_start <= 1;
                    if (n_busy) begin n_start <= 0; state <= S_NORM_BC_W; end
                end

                S_NORM_BC_W: if (n_valid) begin
                    {uBCx, uBCy, uBCz} <= {n_nx, n_ny, n_nz};
                    inv_len_bc <= n_inv_len; // Latch 1/|BC|
                    state <= S_ACOS;
                end

                S_ACOS: begin
                    a_in <= qmult(uBAx, uBCx) + qmult(uBAy, uBCy) + qmult(uBAz, uBCz);
                    a_start <= 1;
                    if (a_busy) begin a_start <= 0; state <= S_ACOS_W; end
                end

                // S_ACOS_W: if (a_valid) begin
                //     d_theta <= theta0 - a_out;
                //     state <= S_CROSS_1;
                //     // The display in the NEXT state will now show the correct d_theta
                // end
                S_ACOS_W: if (a_valid) begin
                    d_theta <= theta0 - a_out;
                    
                    // --- DEBUG LOG: THE HEART OF THE PROBLEM ---
                    $display("[%0t] [ANGLE_CORE] Target: %f | Actual: %f | Delta: %f | K: %f", 
                             $time, 
                             $itor(theta0)/65536.0, 
                             $itor(a_out)/65536.0, 
                             $itor(theta0 - a_out)/65536.0,
                             $itor(k_theta)/65536.0);
                             
                    state <= S_CROSS_1;
                end

                // --- REPAIRED CROSS LOGIC ---
                S_CROSS_1: begin
                    // Find Plane Normal: Normal = BA x BC
                    cp_v1x <= uBAx; cp_v1y <= uBAy; cp_v1z <= uBAz;
                    cp_v2x <= uBCx; cp_v2y <= uBCy; cp_v2z <= uBCz;
                    state <= S_CROSS_A;
                    // $display("[%0t] [CORE] CROSS_1: Inputs BA/BC set for plane normal.", $time);
                end

                S_CROSS_A: begin
                    // 1. CHECK FOR COLLINEARITY (Straight Line)
                    if (cp_rx == 32'sd0 && cp_ry == 32'sd0 && cp_rz == 32'sd0) begin
                        // If atoms are straight, we must "invent" a normal vector.
                        // We test if the molecule is aligned with Z. If not, use Z.
                        // If it IS aligned with Z, use Y.
                        if (uBAx == 32'sd0 && uBAy == 32'sd0) begin
                            {nx, ny, nz} <= {32'sd0, 32'h00010000, 32'sd0}; // Kick Y
                        end else begin
                            {nx, ny, nz} <= {32'sd0, 32'sd0, 32'h00010000}; // Kick Z
                        end
                    end else begin
                        // Standard Case: The atoms already define a clear plane.
                        {nx, ny, nz} <= {cp_rx, cp_ry, cp_rz};
                    end
                    
                    // 2. SET INPUTS FOR FORCE A DIRECTION (Normal x BA)
                    // We use the newly determined nx/ny/nz to ensure cp_v1 isn't zero.
                    cp_v1x <= (cp_rx == 32'sd0 && cp_ry == 32'sd0 && cp_rz == 32'sd0) ? 
                              ((uBAx == 32'sd0 && uBAy == 32'sd0) ? 32'sd0 : 32'sd0) : cp_rx;
                    
                    cp_v1y <= (cp_rx == 32'sd0 && cp_ry == 32'sd0 && cp_rz == 32'sd0) ? 
                              ((uBAx == 32'sd0 && uBAy == 32'sd0) ? 32'h00010000 : 32'sd0) : cp_ry;
                    
                    cp_v1z <= (cp_rx == 32'sd0 && cp_ry == 32'sd0 && cp_rz == 32'sd0) ? 
                              ((uBAx == 32'sd0 && uBAy == 32'sd0) ? 32'sd0 : 32'h00010000) : cp_rz;
                    
                    cp_v2x <= uBAx; 
                    cp_v2y <= uBAy; 
                    cp_v2z <= uBAz;

                    state <= S_CROSS_C;
                end

                S_CROSS_C: begin
                    // Apply 1/r_ba scaling to Force A
                    // fax = k_theta * d_theta * cp_rx * (1/|BA|)
                    fax <= qmult(inv_len_ba, qmult(k_theta, qmult(d_theta, cp_rx)));
                    fay <= qmult(inv_len_ba, qmult(k_theta, qmult(d_theta, cp_ry)));
                    faz <= qmult(inv_len_ba, qmult(k_theta, qmult(d_theta, cp_rz)));

                    cp_v1x <= uBCx; cp_v1y <= uBCy; cp_v1z <= uBCz;
                    cp_v2x <= nx;   cp_v2y <= ny;   cp_v2z <= nz;
                    state <= S_OUT;
                end

                S_OUT: begin
                    // 1. Latch Force C (as we did before)
                    fcx <= qmult(inv_len_bc, qmult(k_theta, qmult(d_theta, cp_rx)));
                    fcy <= qmult(inv_len_bc, qmult(k_theta, qmult(d_theta, cp_ry)));
                    fcz <= qmult(inv_len_bc, qmult(k_theta, qmult(d_theta, cp_rz)));
                    
                    // 2. NEW: Calculate Vertex Force (B) to conserve momentum
                    // Fb = -(Fa + Fc)
                    fbx <= -(fax + qmult(inv_len_bc, qmult(k_theta, qmult(d_theta, cp_rx))));
                    fby <= -(fay + qmult(inv_len_bc, qmult(k_theta, qmult(d_theta, cp_ry))));
                    fbz <= -(faz + qmult(inv_len_bc, qmult(k_theta, qmult(d_theta, cp_rz))));
                    
                    valid_out <= 1; 
                    busy <= 0;
                    state <= S_IDLE;
                end
            endcase
        end
    end
endmodule