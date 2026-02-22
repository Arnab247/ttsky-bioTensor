`default_nettype none

module reciprocal_wrapper (
    input  wire clk,
    input  wire rst_n,
    input  wire valid_in,
    input  wire signed [31:0] x_in,
    output wire signed [31:0] y_out,
    output wire valid_out
);

    // 1. ABSOLUTE VALUE & SIGN
    wire sign = x_in[31];
    wire [31:0] abs_x = sign ? -x_in : x_in;

    // 2. PRIORITY ENCODER (Find Leading 1)
    reg [4:0] lz;
    always @(*) begin
        if      (abs_x[30]) lz = 0;
        else if (abs_x[29]) lz = 1;
        else if (abs_x[28]) lz = 2;
        else if (abs_x[27]) lz = 3;
        else if (abs_x[26]) lz = 4;
        else if (abs_x[25]) lz = 5;
        else if (abs_x[24]) lz = 6;
        else if (abs_x[23]) lz = 7;
        else if (abs_x[22]) lz = 8;
        else if (abs_x[21]) lz = 9;
        else if (abs_x[20]) lz = 10;
        else if (abs_x[19]) lz = 11;
        else if (abs_x[18]) lz = 12;
        else if (abs_x[17]) lz = 13;
        else if (abs_x[16]) lz = 14; 
        else if (abs_x[15]) lz = 15;
        else if (abs_x[14]) lz = 16;
        else if (abs_x[13]) lz = 17;
        else if (abs_x[12]) lz = 18;
        else if (abs_x[11]) lz = 19;
        else if (abs_x[10]) lz = 20;
        else if (abs_x[9])  lz = 21;
        else if (abs_x[8])  lz = 22;
        else if (abs_x[7])  lz = 23;
        else if (abs_x[6])  lz = 24;
        else if (abs_x[5])  lz = 25;
        else if (abs_x[4])  lz = 26;
        else if (abs_x[3])  lz = 27;
        else if (abs_x[2])  lz = 28;
        else if (abs_x[1])  lz = 29;
        else if (abs_x[0])  lz = 30;
        else                lz = 31; // Failsafe for x = 0
    end

    // 3. NORMALIZE INTO [1.0, 2.0)
    // Left-shifting by `lz` puts the leading '1' at bit 30.
    // Right-shifting by 14 moves it to bit 16, which is EXACTLY 1.0 in Q16.16!
    wire [31:0] x_norm = (abs_x << lz) >> 14;

    // 4. CORE MATH (Your rock-solid 6-stage unit)
    wire [31:0] y_norm;
    wire v_out_core;
    
    reciprocal_unit core_math (
        .clk(clk),
        .rst_n(rst_n),
        .valid_in(valid_in),
        .x(x_norm),
        .y_out(y_norm),
        .valid_out(v_out_core)
    );

    // 5. SHIFT PIPELINE
    // The denormalizer needs the shift count to arrive at the exact 
    // same time as the data (6 clock cycles later).
    reg [4:0] lz_pipe [0:5];
    reg       sign_pipe [0:5];
    integer k;
    
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            for (k = 0; k < 6; k = k + 1) begin
                lz_pipe[k]   <= 0;
                sign_pipe[k] <= 0;
            end
        end else begin
            lz_pipe[0]   <= lz;
            sign_pipe[0] <= sign;
            for (k = 1; k < 6; k = k + 1) begin
                lz_pipe[k]   <= lz_pipe[k-1];
                sign_pipe[k] <= sign_pipe[k-1];
            end
        end
    end

    // 6. DENORMALIZE
    // We use a massive 64-bit register here so shifting doesn't destroy bits.
    // If we shifted the input, we apply the EXACT same shift to the output math.
    wire [63:0] y_expanded = {32'b0, y_norm};
    wire [63:0] y_shifted  = (y_expanded << lz_pipe[5]) >> 14;
    
    // Saturation logic (Prevents true mathematically impossible Q16.16 values from wrapping)
    wire [31:0] y_unsigned = (y_shifted > 64'h000000007FFFFFFF) ? 32'h7FFFFFFF : y_shifted[31:0];
    
    assign y_out = sign_pipe[5] ? -y_unsigned : y_unsigned;
    assign valid_out = v_out_core;

endmodule