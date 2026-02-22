`default_nettype none

module acos_poly (
    input  wire        clk,
    input  wire        rst_n,
    input  wire        start,
    input  wire signed [31:0] x_in,
    output reg  signed [31:0] theta_out,
    output reg         valid_out,
    output reg         busy
);

    // -------------------------------------------------------------------------
    // COEFFICIENT LOOKUP TABLES (Fixed-point Q16.16)
    // -------------------------------------------------------------------------
    // SET 1: Center Maclaurin (For -0.75 to +0.75)
    localparam signed [31:0] MAC_C0 = 32'h0001921F; // 1.5707
    localparam signed [31:0] MAC_C1 = 32'hFFFF030A; // -1.0
    localparam signed [31:0] MAC_C2 = 32'h00000000; // 0.0
    localparam signed [31:0] MAC_C3 = 32'hFFFFD555; // -0.1666

    // SET 2: Edge Linear Fits (Slope = -2.09)
    localparam signed [31:0] EDGE_C1     = 32'hFFFDE8F6; // -2.09
    localparam signed [31:0] EDGE_POS_C0 = 32'h000250A3; // 2.315
    localparam signed [31:0] EDGE_NEG_C0 = 32'h0000D374; // 0.826

    // Thresholds
    localparam signed [31:0] UPPER_CLAMP = 32'h0000F000; //  0.9375
    localparam signed [31:0] UPPER_SHOULDER= 32'h0000C000; //  0.75
    localparam signed [31:0] LOWER_SHOULDER= 32'hFFFF4000; // -0.75
    localparam signed [31:0] LOWER_CLAMP = 32'hFFFF1000; // -0.9375

    // Constants
    localparam signed [31:0] PI_RADS = 32'h0003243F; // 3.14159
    localparam signed [31:0] ZERO    = 32'h00000000;

    reg [2:0] stage;
    reg signed [31:0] x_reg, acc;
    reg override_active;
    reg signed [31:0] override_val;

    // Dynamically selected coefficients
    reg signed [31:0] curr_c0, curr_c1, curr_c2, curr_c3;

    function automatic signed [31:0] qmult(input signed [31:0] a, input signed [31:0] b);
        reg signed [63:0] temp;
        begin
            temp = $signed({{32{a[31]}}, a}) * $signed({{32{b[31]}}, b});
            qmult = temp[47:16];
        end
    endfunction

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            stage <= 0; busy <= 0; valid_out <= 0; theta_out <= 0;
            override_active <= 0; override_val <= 0; x_reg <= 0; acc <= 0;
            curr_c0 <= 0; curr_c1 <= 0; curr_c2 <= 0; curr_c3 <= 0;
        end else begin
            case (stage)
                0: begin
                    valid_out <= 0;
                    if (start) begin
                        busy <= 1;
                        x_reg <= x_in;
                        override_active <= 0;

                        // --- ZERO-CYCLE HARDWARE ROUTER ---
                        if (x_in >= UPPER_CLAMP) begin
                            override_active <= 1; override_val <= ZERO;
                        end 
                        else if (x_in <= LOWER_CLAMP) begin
                            override_active <= 1; override_val <= PI_RADS;
                        end 
                        else if (x_in > UPPER_SHOULDER) begin
                            // Positive Shoulder (Linear Fit)
                            curr_c0 <= EDGE_POS_C0; curr_c1 <= EDGE_C1; 
                            curr_c2 <= ZERO; curr_c3 <= ZERO;
                            acc <= ZERO; // Start with 0 since C3 is 0
                        end
                        else if (x_in < LOWER_SHOULDER) begin
                            // Negative Shoulder (Linear Fit)
                            curr_c0 <= EDGE_NEG_C0; curr_c1 <= EDGE_C1; 
                            curr_c2 <= ZERO; curr_c3 <= ZERO;
                            acc <= ZERO; 
                        end
                        else begin
                            // Center Zone (Maclaurin Curve)
                            curr_c0 <= MAC_C0; curr_c1 <= MAC_C1; 
                            curr_c2 <= MAC_C2; curr_c3 <= MAC_C3;
                            acc <= MAC_C3;
                        end
                        
                        stage <= 1;
                    end else busy <= 0;
                end
                1: begin acc <= curr_c2 + qmult(x_reg, acc); stage <= 2; end
                2: begin acc <= curr_c1 + qmult(x_reg, acc); stage <= 3; end
                3: begin theta_out <= curr_c0 + qmult(x_reg, acc); stage <= 4; end
                4: begin
                    if (override_active) theta_out <= override_val;
                    valid_out <= 1; busy <= 0; stage <= 0;
                end
                default: stage <= 0;
            endcase
        end
    end
endmodule