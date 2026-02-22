`default_nettype none

module vector_normalizer (
    input  wire        clk,
    input  wire        rst_n,
    input  wire        start,
    
    input  wire signed [31:0] vx, vy, vz,
    
    output reg  signed [31:0] nx, ny, nz,
    output reg         valid_out,
    output reg signed [31:0] inv_mag, // Added output for 1/|V| in Q16.16
    output reg         busy
);

    // FSM States
    localparam S_IDLE  = 0;
    localparam S_MAC   = 1;
    localparam S_SQRT  = 2;
    localparam S_DIV   = 3;
    localparam S_NORM  = 4;

    reg [2:0] state;
    reg [63:0] sq_sum;
    reg [63:0] res;
    reg [63:0] curr_bit;
    // reg signed [31:0] inv_mag;

    // ==========================================
    // EXPLICIT COMBINATIONAL 64-BIT CASTING
    // (Prevents simulation race conditions)
    // ==========================================
    wire signed [63:0] w_vx = {{32{vx[31]}}, vx};
    wire signed [63:0] w_vy = {{32{vy[31]}}, vy};
    wire signed [63:0] w_vz = {{32{vz[31]}}, vz};

    wire signed [63:0] sq_x = w_vx * w_vx;
    wire signed [63:0] sq_y = w_vy * w_vy;
    wire signed [63:0] sq_z = w_vz * w_vz;
    // ==========================================

    function automatic signed [31:0] qmult(input signed [31:0] a, input signed [31:0] b);
        reg signed [63:0] temp;
        begin
            temp = $signed({{32{a[31]}}, a}) * $signed({{32{b[31]}}, b});
            qmult = temp[47:16];
        end
    endfunction

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            state <= S_IDLE;
            valid_out <= 0; busy <= 0;
            nx <= 0; ny <= 0; nz <= 0;
        end else begin
            valid_out <= 0; // Default off
            
            case (state)
                S_IDLE: begin
                    if (start) begin
                        busy <= 1;
                        state <= S_MAC;
                    end
                end
                
                S_MAC: begin
                    // Read the stable combinational wires
                    sq_sum <= sq_x + sq_y + sq_z;
                    res <= 0;
                    curr_bit <= 64'h4000000000000000;
                    state <= S_SQRT;
                end
                
                S_SQRT: begin
                    // Hardware Digit-by-Digit Square Root
                    if (curr_bit > 0) begin
                        if (sq_sum >= res + curr_bit) begin
                            sq_sum <= sq_sum - (res + curr_bit);
                            res <= (res >> 1) + curr_bit;
                        end else begin
                            res <= res >> 1;
                        end
                        curr_bit <= curr_bit >> 2;
                    end else begin
                        state <= S_DIV;
                    end
                end
                
                S_DIV: begin
                    // Hardware Division (Yields exact Q16.16 scalar)
                    if (res == 0) begin
                        inv_mag <= 0;
                    end else begin
                        inv_mag <= (64'h0000000100000000) / res; 
                    end
                    state <= S_NORM;
                end
                
                S_NORM: begin
                    nx <= qmult(vx, inv_mag);
                    ny <= qmult(vy, inv_mag);
                    nz <= qmult(vz, inv_mag);
                    valid_out <= 1;
                    state <= S_IDLE;
                    busy <= 0;
                end
            endcase
        end
    end
endmodule