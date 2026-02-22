`default_nettype none

module bonded_force_core (
    input  wire        clk,
    input  wire        rst_n,
    input  wire        start,
    
    // Atom Pair Positions
    input  wire signed [31:0] x1, y1, z1,
    input  wire signed [31:0] x2, y2, z2,
    
    // Bond Parameters (Q16.16)
    input  wire signed [31:0] r0,  // Equilibrium length
    input  wire signed [31:0] k,   // Spring constant
    
    // Output Forces (Applied to Atom 1. Atom 2 gets the negative of this)
    output reg signed [31:0] fx, fy, fz,
    output reg         valid_out,
    output reg         busy
);

    localparam S_IDLE  = 0;
    localparam S_DELTA = 1;
    localparam S_SQRT  = 2;
    localparam S_DIV   = 3;
    localparam S_FORCE = 4;

    reg [2:0] state;
    reg signed [31:0] dx, dy, dz;
    reg [63:0] sq_sum;
    reg [63:0] r_res;
    reg [63:0] curr_bit;
    
    reg signed [31:0] r;        // Actual distance
    reg signed [31:0] inv_r;    // 1.0 / distance
    reg signed [31:0] f_scalar; // Magnitude of force

    // Explicit 64-bit casts for safe multiplication
    reg signed [63:0] dx64, dy64, dz64;

    function signed [31:0] qmult(input signed [31:0] a, input signed [31:0] b);
        reg signed [63:0] temp;
        begin
            temp = $signed({{32{a[31]}}, a}) * $signed({{32{b[31]}}, b});
            qmult = temp[47:16];
        end
    endfunction

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            state <= S_IDLE;
            busy <= 0; valid_out <= 0;
            fx <= 0; fy <= 0; fz <= 0;
        end else begin
            valid_out <= 0;
            
            case (state)
                S_IDLE: begin
                    if (start) begin
                        busy <= 1;
                        state <= S_DELTA;
                    end
                end

                S_DELTA: begin
                    dx = x2 - x1;
                    dy = y2 - y1;
                    dz = z2 - z1;
                    
                    dx64 = $signed(dx); 
                    dy64 = $signed(dy); 
                    dz64 = $signed(dz);
                    
                    sq_sum <= (dx64 * dx64) + (dy64 * dy64) + (dz64 * dz64);
                    r_res <= 0;
                    curr_bit <= 64'h4000000000000000;
                    state <= S_SQRT;
                end

                S_SQRT: begin
                    // Hardware Square Root to find distance 'r'
                    if (curr_bit > 0) begin
                        if (sq_sum >= r_res + curr_bit) begin
                            sq_sum <= sq_sum - (r_res + curr_bit);
                            r_res <= (r_res >> 1) + curr_bit;
                        end else begin
                            r_res <= r_res >> 1;
                        end
                        curr_bit <= curr_bit >> 2;
                    end else begin
                        r <= r_res[31:0]; // Save actual distance
                        state <= S_DIV;
                    end
                end

                S_DIV: begin
                    if (r == 0) begin
                        inv_r <= 0;
                    end else begin
                        inv_r <= (64'h0000000100000000) / r_res; 
                    end
                    state <= S_FORCE;
                end

                S_FORCE: begin
                    // F_scalar = 2 * k * (r - r0)
                    // If r > r0 (stretched), F_scalar is positive (attractive)
                    f_scalar = qmult(k <<< 1, r - r0);
                    
                    // Project onto axes: Fx = F_scalar * (dx / r)
                    fx <= qmult(f_scalar, qmult(dx, inv_r));
                    fy <= qmult(f_scalar, qmult(dy, inv_r));
                    fz <= qmult(f_scalar, qmult(dz, inv_r));
                    
                    valid_out <= 1;
                    busy <= 0;
                    state <= S_IDLE;
                end
            endcase
        end
    end
endmodule