`default_nettype none

module cordic_atan2 (
    input  wire        clk,
    input  wire        rst_n,
    input  wire        start,
    input  wire signed [31:0] x_in,
    input  wire signed [31:0] y_in,
    
    output reg signed [31:0] theta_out,
    output reg         valid_out,
    output reg         busy
);

    // Q16.16 Pi Constants
    localparam signed [31:0] PI_OVER_2  = 32'h0001921F; 
    localparam signed [31:0] N_PI_OVER_2= 32'hFFFE6DE1; 

    // LUT
    reg signed [31:0] atan_lut [0:15];
    initial begin
        atan_lut[0] = 32'h0000C90F; atan_lut[1] = 32'h000076B1; atan_lut[2] = 32'h00003EB6;
        atan_lut[3] = 32'h00001FD5; atan_lut[4] = 32'h00000FFA; atan_lut[5] = 32'h000007FF;
        atan_lut[6] = 32'h00000400; atan_lut[7] = 32'h00000200; atan_lut[8] = 32'h00000100;
        atan_lut[9] = 32'h00000080; atan_lut[10]= 32'h00000040; atan_lut[11]= 32'h00000020;
        atan_lut[12]= 32'h00000010; atan_lut[13]= 32'h00000008; atan_lut[14]= 32'h00000004;
        atan_lut[15]= 32'h00000002;
    end

    reg [4:0]  step;
    reg signed [31:0] x_reg, y_reg, z_reg;
    
    // FIX 1: Move shifters to combinational logic (WIRE)
    // This guarantees they update instantly when x_reg/y_reg change
    wire signed [31:0] x_shift = x_reg >>> step;
    wire signed [31:0] y_shift = y_reg >>> step;

    localparam S_IDLE = 0, S_PRE = 1, S_CALC = 2, S_DONE = 3;
    reg [1:0] state;

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            state <= S_IDLE;
            valid_out <= 0; busy <= 0;
            theta_out <= 0; step <= 0;
            x_reg <= 0; y_reg <= 0; z_reg <= 0;
        end else begin
            case (state)
                S_IDLE: begin
                    valid_out <= 0;
                    if (start) begin
                        busy <= 1;
                        step <= 0;
                        state <= S_PRE; // Move to explicit pre-rotation state
                    end
                end

                S_PRE: begin
                    // FIX 2: Explicit Pre-Rotation State
                    // This handles the +/- 90 degree rotation cleanly before iterating
                    if (x_in >= 0) begin
                        x_reg <= x_in; 
                        y_reg <= y_in; 
                        z_reg <= 0;
                    end else if (y_in >= 0) begin // Quadrant II
                        x_reg <= y_in; 
                        y_reg <= -x_in; 
                        z_reg <= PI_OVER_2;
                    end else begin                // Quadrant III
                        x_reg <= -y_in; 
                        y_reg <= x_in; 
                        z_reg <= N_PI_OVER_2;
                    end
                    state <= S_CALC;
                end

                S_CALC: begin
                    // Vectoring Mode: Drive Y to zero
                    // Note: We use the 'wire' x_shift/y_shift here
                    if (y_reg[31] == 0) begin // If Y is positive
                        x_reg <= x_reg + y_shift;
                        y_reg <= y_reg - x_shift;
                        z_reg <= z_reg + atan_lut[step];
                    end else begin            // If Y is negative
                        x_reg <= x_reg - y_shift;
                        y_reg <= y_reg + x_shift;
                        z_reg <= z_reg - atan_lut[step];
                    end

                    if (step == 15) begin
                        state <= S_DONE;
                    end else begin
                        step <= step + 1;
                    end
                end

                S_DONE: begin
                    // FIX 3: Clean handshake
                    theta_out <= z_reg;
                    valid_out <= 1;
                    busy <= 0;
                    state <= S_IDLE; // Go to IDLE immediately
                end
            endcase
        end
    end
endmodule