`default_nettype none

module md_system_top (
    input  wire        clk,
    input  wire        rst_n,
    input  wire        start_run,
    input  wire [15:0] max_iters,
    input  wire [5:0]  num_atoms,

    input  wire        load_en,
    input  wire [5:0]  load_addr,
    input  wire signed [31:0] load_x, load_y, load_z,
    input  wire [4:0]  load_res_id,
    input  wire [3:0]  load_atom_idx,

    output reg         done,
    output reg  [15:0] iter_count
);

    // =========================================================================
    // STEP SIZE (Learning Rate / Alpha = 0.001)
    // =========================================================================
    // localparam signed [31:0] STEP_SIZE = 32'h00000041;
    // REPLACE: localparam signed [31:0] STEP_SIZE = 32'h00000041;
    // WITH:
    reg signed [31:0] current_step_size;
    localparam signed [31:0] INITIAL_STEP  = 32'h00000200; // ~0.01 (10x larger)
    localparam signed [31:0] COOLING_FACTOR = 32'h0000F000; // ~0.80 in Q16.16

    // =========================================================================
    // FSM STATE ENCODING (Two-Phase Architecture)
    // =========================================================================
    localparam S_IDLE         = 0,
               S_ITER_START   = 1,
               
               // PHASE 1: Bonded Sliding Window
               S_BND_FETCH    = 2,
               S_BND_LOOKUP   = 3,
               S_BND_EVAL     = 4,
               S_BND_WAIT     = 5,
               S_BND_ACCUM    = 6,
               S_BND_NEXT     = 7,
               
               // PHASE 2: Non-Bonded Cloud (Nested Loop)
               S_NB_START     = 8,
               S_NB_INNER     = 9,
               S_NB_FETCH     = 10,
               S_NB_LOOKUP    = 11,
               S_NB_FEED      = 12,
               S_NB_DRAIN     = 13,
               
               // PHASE 3: Gradient Update
               S_APPLY_FETCH  = 14,
               S_APPLY_UPDATE = 15,
               S_APPLY_NEXT   = 16;

    reg [4:0] state;
    reg [5:0] scan_idx;
    reg [5:0] apply_idx;

    // Non-Bonded loop counters
    reg [5:0] nb_i, nb_j;
    reg [15:0] nb_inflight; // Tracks outstanding pipeline calculations

    // =========================================================================
    // FORCE ACCUMULATOR ARRAYS
    // =========================================================================
    reg signed [31:0] acc_fx [0:63];
    reg signed [31:0] acc_fy [0:63];
    reg signed [31:0] acc_fz [0:63];
    integer i;

    // =========================================================================
    // ATOM REGISTER FILE WIRING
    // =========================================================================
    reg        we;
    reg  [5:0] w_addr;
    reg  signed [31:0] w_x, w_y, w_z;
    reg  [5:0] r_a, r_b, r_c, r_d;
    wire signed [31:0] xa, ya, za, xb, yb, zb, xc, yc, zc, xd, yd, zd;
    wire [4:0] ra_res_id;
    wire [3:0] ra_atom_idx;

    atom_regfile u_memory (
        .clk      (clk),
        .rst_n    (rst_n),
        .we       (load_en | we),
        .w_addr   (load_en ? load_addr   : w_addr),
        .w_x      (load_en ? load_x      : w_x),
        .w_y      (load_en ? load_y      : w_y),
        .w_z      (load_en ? load_z      : w_z),
        .w_res_id (load_res_id),
        .w_atom_idx(load_atom_idx),
        .r_addr_a (r_a), .r_addr_b(r_b), .r_addr_c(r_c), .r_addr_d(r_d),
        .xa(xa), .ya(ya), .za(za),
        .xb(xb), .yb(yb), .zb(zb),
        .xc(xc), .yc(yc), .zc(zc),
        .xd(xd), .yd(yd), .zd(zd),
        .ra_res_id  (ra_res_id),
        .ra_atom_idx(ra_atom_idx)
    );

    // =========================================================================
    // PARAMETER RAM & DATABASES (From Teammate)
    // =========================================================================
    wire [31:0] db_param_out;
    wire [15:0] current_type_id = db_param_out[15:0];
    wire [31:0] sigma, epsilon;
    wire signed [31:0] mem_r0, mem_kb, mem_theta0, mem_k_theta, mem_phi0, mem_k_phi;
    wire [3:0]  mem_n_period;
    wire signed [31:0] mem_q_a, mem_q_d;

    residue_database u_rdb (.residue_id(ra_res_id), .atom_index(ra_atom_idx), .param_out(db_param_out));
    atom_type_table u_att (.type_id(current_type_id), .sigma_q16(sigma), .epsilon_q16(epsilon));
    parameter_ram u_param_memory (
        .clk(clk), 
        // Multiplex RAM address based on phase (bonded vs non-bonded)
        .addr({7'b0, scan_idx[2:0]}), 
        .r0_out(mem_r0),
        .kb_out(mem_kb),
        .theta0_out(mem_theta0),
        .k_theta_out(mem_k_theta),
        .phi0_out(mem_phi0), .k_phi_out(mem_k_phi), .n_period_out(mem_n_period),
        .q_a_out(mem_q_a), .q_d_out(mem_q_d)
    );

    // Active Latch Registers
    reg signed [31:0] active_r0, active_kb, active_theta0, active_k_theta;
    reg signed [31:0] active_phi0, active_k_phi, active_q_a, active_q_d;
    reg signed [31:0] active_sigma_sq, active_eps_x24;
    reg [3:0] active_n_period;

    // =========================================================================
    // PHYSICS CORES
    // =========================================================================
    reg bnd_start, bnd_done; wire bnd_valid, bnd_busy;
    wire signed [31:0] f_bnd_ax, f_bnd_ay, f_bnd_az;
    wire signed [31:0] f_bnd_bx = -f_bnd_ax, f_bnd_by = -f_bnd_ay, f_bnd_bz = -f_bnd_az;
    bonded_force_core u_bond (.clk(clk), .rst_n(rst_n), .start(bnd_start), .x1(xa), .y1(ya), .z1(za), .x2(xb), .y2(yb), .z2(zb), .r0(active_r0), .k(active_kb), .fx(f_bnd_ax), .fy(f_bnd_ay), .fz(f_bnd_az), .valid_out(bnd_valid), .busy(bnd_busy));

    reg ang_start, ang_done; wire ang_valid, ang_busy;
    wire signed [31:0] f_ang_ax, f_ang_ay, f_ang_az, f_ang_cx, f_ang_cy, f_ang_cz;

    // --- Angle Force Output Wires ---
    wire signed [31:0] fax, fay, faz;
    wire signed [31:0] fbx, fby, fbz;
    wire signed [31:0] fcx, fcy, fcz;


    // angle_force_core u_angle (.clk(clk), .rst_n(rst_n), .start(ang_start), .xa(xa), .ya(ya), .za(za), .xb(xb), .yb(yb), .zb(zb), .xc(xc), .yc(yc), .zc(zc), .theta0(active_theta0), .k_theta(active_k_theta), .fax(f_ang_ax), .fay(f_ang_ay), .faz(f_ang_az), .fcx(f_ang_cx), .fcy(f_ang_cy), .fcz(f_ang_cz), .valid_out(ang_valid), .busy(ang_busy));
    angle_force_core u_angle (
        .clk(clk),
        .rst_n(rst_n),
        .start(ang_start),
        .xa(xa), .ya(ya), .za(za),
        .xb(xb), .yb(yb), .zb(zb),
        .xc(xc), .yc(yc), .zc(zc),
        .theta0(active_theta0),
        .k_theta(active_k_theta),
        .fax(fax), .fay(fay), .faz(faz), // These now bind to the wires above
        .fbx(fbx), .fby(fby), .fbz(fbz),
        .fcx(fcx), .fcy(fcy), .fcz(fcz),
        .valid_out(ang_valid),
        .busy(ang_busy)
    );




    reg dih_start, dih_done; wire dih_valid, dih_busy;
    // wire signed [31:0] f_dih_ax, f_dih_ay, f_dih_az, f_dih_bx, f_dih_by, f_dih_bz, f_dih_cx, f_dih_cy, f_dih_cz, f_dih_dx, f_dih_dy, f_dih_dz;
    wire signed [31:0] f_dih_ax, f_dih_ay, f_dih_az;
    wire signed [31:0] f_dih_bx, f_dih_by, f_dih_bz;
    wire signed [31:0] f_dih_cx, f_dih_cy, f_dih_cz;
    wire signed [31:0] f_dih_dx, f_dih_dy, f_dih_dz;
    dihedral_force_core u_dihedral (.clk(clk), .rst_n(rst_n), .start(dih_start), .xa(xa), .ya(ya), .za(za), .xb(xb), .yb(yb), .zb(zb), .xc(xc), .yc(yc), .zc(zc), .xd(xd), .yd(yd), .zd(zd), .phi0(active_phi0), .k_phi(active_k_phi), .n_period(active_n_period), .fax(f_dih_ax), .fay(f_dih_ay), .faz(f_dih_az), .fbx(f_dih_bx), .fby(f_dih_by), .fbz(f_dih_bz), .fcx(f_dih_cx), .fcy(f_dih_cy), .fcz(f_dih_cz), .fdx(f_dih_dx), .fdy(f_dih_dy), .fdz(f_dih_dz), .valid_out(dih_valid), .busy(dih_busy));

    reg nb_valid_in;
    wire nb_valid_out;
    wire signed [31:0] f_nb_x, f_nb_y, f_nb_z;
    non_bonded_pipeline u_nbpipe (.clk(clk), .rst_n(rst_n), .valid_in(nb_valid_in), .xi(xa), .yi(ya), .zi(za), .xj(xd), .yj(yd), .zj(zd), .q_i(active_q_a), .q_j(active_q_d), .sigma_sq(active_sigma_sq), .eps_x24(active_eps_x24), .valid_out(nb_valid_out), .fx_nb(f_nb_x), .fy_nb(f_nb_y), .fz_nb(f_nb_z));




    // 16-bit LFSR for pseudo-random noise
    reg [15:0] lfsr;
    wire signed [31:0] random_force;
    
    // Feedback polynomial for a maximal-period 16-bit LFSR
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) lfsr <= 16'hACE1; // Seed value
        else lfsr <= {lfsr[14:0], lfsr[15] ^ lfsr[13] ^ lfsr[12] ^ lfsr[10]};
    end

    // Scale the 16-bit random number to a small Q16.16 signed force
    // This creates a "kick" between -0.5 and +0.5
    assign random_force = {{16{lfsr[15]}}, lfsr};



    // =========================================================================
    // PIPELINE ALIGNMENT SHIFT REGISTERS (For Newton's 3rd Law in Phase 2)
    // =========================================================================
    reg [5:0] nb_i_sr [0:17];
    reg [5:0] nb_j_sr [0:17];
    always @(posedge clk) begin
        nb_i_sr[0] <= nb_i;
        nb_j_sr[0] <= nb_j;
        for (integer k=1; k<=17; k=k+1) begin
            nb_i_sr[k] <= nb_i_sr[k-1];
            nb_j_sr[k] <= nb_j_sr[k-1];
        end
    end

    // =========================================================================
    // MATH HELPERS
    // =========================================================================
    function automatic signed [31:0] qmult(input signed [31:0] a, input signed [31:0] b);
        reg signed [63:0] temp; begin temp = $signed({{32{a[31]}}, a}) * $signed({{32{b[31]}}, b}); qmult = temp[47:16]; end
    endfunction

    function automatic signed [31:0] clamp_move(input signed [31:0] move);
        // 0x00008000 = 0.5 in Q16.16
        localparam signed [31:0] LIMIT = 32'h00008000;
        begin
            if (move > LIMIT)  clamp_move = LIMIT;
            else if (move < -LIMIT) clamp_move = -LIMIT;
            else clamp_move = move;
        end
    endfunction
    

    // =========================================================================
    // MAIN MULTI-PHASE FSM
    // =========================================================================
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            state <= S_IDLE; done <= 0; iter_count <= 0; scan_idx <= 0; apply_idx <= 0;
            we <= 0; w_addr <= 0; w_x <= 0; w_y <= 0; w_z <= 0; r_a <= 0; r_b <= 0; r_c <= 0; r_d <= 0;
            {bnd_start, ang_start, dih_start} <= 3'b0; {bnd_done, ang_done, dih_done} <= 3'b0; nb_valid_in <= 0; nb_inflight <= 0;
            for (i = 0; i < 64; i = i + 1) begin acc_fx[i] <= 0; acc_fy[i] <= 0; acc_fz[i] <= 0; end
        end else begin
            
            // INDEPENDENT ACCUMULATION: Always listen for the NB pipeline finish line
            if (nb_valid_out) begin
                acc_fx[nb_i_sr[17]] <= acc_fx[nb_i_sr[17]] + f_nb_x;
                acc_fy[nb_i_sr[17]] <= acc_fy[nb_i_sr[17]] + f_nb_y;
                acc_fz[nb_i_sr[17]] <= acc_fz[nb_i_sr[17]] + f_nb_z;
                // Newton's 3rd Law: Equal and opposite force on Atom J
                acc_fx[nb_j_sr[17]] <= acc_fx[nb_j_sr[17]] - f_nb_x;
                acc_fy[nb_j_sr[17]] <= acc_fy[nb_j_sr[17]] - f_nb_y;
                acc_fz[nb_j_sr[17]] <= acc_fz[nb_j_sr[17]] - f_nb_z;
                nb_inflight <= nb_inflight - 16'd1;
            end

            case (state)
                S_IDLE: begin
                    done <= 0;
                    if (start_run) begin iter_count <= 0; state <= S_ITER_START; end
                end

                S_ITER_START: begin
                    for (i = 0; i < 64; i = i + 1) begin 
                        acc_fx[i] <= 0; acc_fy[i] <= 0; acc_fz[i] <= 0; 
                    end
                    scan_idx <= 0;
                    // Keep the annealing step size logic
                    current_step_size <= (iter_count == 0) ? INITIAL_STEP : current_step_size; 
                    state <= S_BND_FETCH;
                end

                // --- PHASE 1: BONDED SLIDING WINDOW ---
                S_BND_FETCH: begin
                    r_a <= scan_idx; r_b <= scan_idx + 6'd1; r_c <= scan_idx + 6'd2; r_d <= scan_idx + 6'd3;
                    bnd_done <= 0; ang_done <= 0; dih_done <= 0; state <= S_BND_LOOKUP;
                end
                S_BND_LOOKUP: begin
                    active_r0 <= mem_r0; active_kb <= mem_kb; active_theta0 <= mem_theta0; active_k_theta <= mem_k_theta;
                    active_phi0 <= mem_phi0; active_k_phi <= mem_k_phi; active_n_period <= mem_n_period;
                    state <= S_BND_EVAL;
                end
                S_BND_EVAL: begin bnd_start <= 1; ang_start <= 1; dih_start <= 1; state <= S_BND_WAIT; end
                S_BND_WAIT: begin
                    bnd_start <= 0; ang_start <= 0; dih_start <= 0;
                    if (bnd_valid) bnd_done <= 1;
                    if (ang_valid) ang_done <= 1;
                    if (dih_valid) dih_done <= 1;
                    if ((bnd_done | bnd_valid) && (ang_done | ang_valid) && (dih_done | dih_valid)) state <= S_BND_ACCUM;
                end
                S_BND_ACCUM: begin
                    // ATOM A: Bond + Angle + Dihedral
                    acc_fx[r_a] <= acc_fx[r_a] + f_bnd_ax + (fax) + f_dih_ax;
                    acc_fy[r_a] <= acc_fy[r_a] + f_bnd_ay + (fay) + f_dih_ay;
                    acc_fz[r_a] <= acc_fz[r_a] + f_bnd_az + (faz) + f_dih_az;

                    // ATOM B: Opposite Bond + Pivot Angle + Dihedral
                    // Note: f_bnd_bx is -f_bnd_ax
                    acc_fx[r_b] <= acc_fx[r_b] + f_bnd_bx + (fbx) + f_dih_bx;
                    acc_fy[r_b] <= acc_fy[r_b] + f_bnd_by + (fby) + f_dih_by;
                    acc_fz[r_b] <= acc_fz[r_b] + f_bnd_bz + (fbz) + f_dih_bz;

                    // ATOM C: Angle + Dihedral
                    acc_fx[r_c] <= acc_fx[r_c] + (fcx) + f_dih_cx;
                    acc_fy[r_c] <= acc_fy[r_c] + (fcy) + f_dih_cy;
                    acc_fz[r_c] <= acc_fz[r_c] + (fcz) + f_dih_cz;

                    // SYMMETRY BREAKER: If Z is perfectly flat, give a tiny push
                    // We check za, zb, zc directly from the regfile wires
                    if (za == 0 && zb == 0 && zc == 0) begin
                        acc_fz[r_b] <= acc_fz[r_b] + (random_force >>> 10);
                    end

                    state <= S_BND_NEXT;
                end
                S_BND_NEXT: begin
                    if (scan_idx + 6'd4 < {1'b0, num_atoms}) begin
                        scan_idx <= scan_idx + 6'd1; state <= S_BND_FETCH;
                    end else begin
                        state <= S_NB_START; // Move to Phase 2!
                    end
                end

                // --- PHASE 2: NON-BONDED O(N^2) CLOUD ---
                S_NB_START: begin nb_i <= 0; nb_j <= 3; nb_inflight <= 0; state <= S_NB_INNER; end
                S_NB_INNER: begin
                    nb_valid_in <= 0; // default off
                    if (nb_i >= num_atoms - 6'd1) begin
                        state <= S_NB_DRAIN; // All pairs fired, wait for pipeline
                    end else if (nb_j >= num_atoms) begin
                        nb_i <= nb_i + 6'd1;
                        nb_j <= nb_i + 6'd4; // Enforce i+3 exclusion jump!
                    end else begin
                        r_a <= nb_i; r_d <= nb_j; state <= S_NB_FETCH;
                    end
                end
                S_NB_FETCH: begin state <= S_NB_LOOKUP; end
                S_NB_LOOKUP: begin
                    active_q_a <= mem_q_a; active_q_d <= mem_q_d; // Fetch charges
                    active_sigma_sq <= sigma; active_eps_x24 <= epsilon;
                    state <= S_NB_FEED;
                end
                S_NB_FEED: begin
                    nb_valid_in <= 1; // Fire the pipeline!
                    nb_inflight <= nb_inflight + 16'd1;
                    nb_j <= nb_j + 6'd1;
                    state <= S_NB_INNER;
                end
                S_NB_DRAIN: begin
                    if (nb_inflight == 0 && !nb_valid_out) begin
                        apply_idx <= 0; state <= S_APPLY_FETCH; // Move to Phase 3!
                    end
                end

                // --- PHASE 3: GRADIENT DESCENT UPDATE ---
                S_APPLY_FETCH: begin r_a <= apply_idx; state <= S_APPLY_UPDATE; end
                
                S_APPLY_UPDATE: begin
                    we     <= 1;
                    w_addr <= apply_idx;
                    
                    // Limit the maximum move per iteration to 0.5 Angstroms (0x00008000)
                    // This prevents a single bad calculation from teleporting an atom to infinity
                    w_x <= xa + clamp_move(qmult(acc_fx[apply_idx], current_step_size));
                    w_y <= ya + clamp_move(qmult(acc_fy[apply_idx], current_step_size));
                    w_z <= za + clamp_move(qmult(acc_fz[apply_idx], current_step_size));
                    
                    state  <= S_APPLY_NEXT;
                end

                S_APPLY_NEXT: begin
                    we <= 0;
                    if (apply_idx < num_atoms - 6'd1) begin
                        apply_idx <= apply_idx + 6'd1; 
                        state <= S_APPLY_FETCH;
                    end else begin
                        iter_count <= iter_count + 16'd1;
                        // COOL THE SYSTEM: current_step_size = current_step_size * 0.8
                        current_step_size <= qmult(current_step_size, COOLING_FACTOR);
                        
                        if (iter_count + 16'd1 >= max_iters) begin 
                            done <= 1; state <= S_IDLE; 
                        end else begin 
                            state <= S_ITER_START; 
                        end
                    end
                end
            endcase
        end
    end

endmodule