`default_nettype none

module reciprocal_unit (
    input  wire clk,
    input  wire rst_n,
    input  wire valid_in,
    input  wire signed [31:0] x,              
    output reg  signed [31:0] y_out,
    output reg  valid_out
);

    // ---------------------------------------------------------
    // MATH HELPER: Q16.16 Multiplier with Convergent Rounding
    // ---------------------------------------------------------
    function automatic signed [31:0] qmult(input signed [31:0] a, input signed [31:0] b);
        reg signed [63:0] temp;
        begin
            temp = $signed({{32{a[31]}}, a}) * $signed({{32{b[31]}}, b});
            // Add 0.5 LSB (offset 0x8000) before truncation for better precision
            temp = temp + 64'h8000;
            qmult = temp[47:16]; 
        end
    endfunction

    localparam signed [31:0] TWO = 32'h00020000;

    wire signed [31:0] y0;
    norm_seed_lut u_seed (.d_in(x), .seed_out(y0)); 

    // --- PIPELINE REGISTERS ---
    reg signed [31:0] y1, y2, y3, y4, y5;
    reg signed [31:0] x_d1, x_d2, x_d3, x_d4, x_d5;
    reg v1, v2, v3, v4, v5; 

    // STAGE 1: Capture x and y0 on the cycle when valid_in is high
    wire signed [31:0] term1 = TWO - qmult(x, y0);
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin y1 <= 0; x_d1 <= 0; v1 <= 0; end 
        else if (valid_in) begin 
            y1 <= qmult(y0, term1); 
            x_d1 <= x; 
            v1 <= 1; 
        end else begin
            v1 <= 0;
        end
    end

    // STAGE 2
    wire signed [31:0] term2 = TWO - qmult(x_d1, y1);
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin y2 <= 0; x_d2 <= 0; v2 <= 0; end 
        else begin y2 <= qmult(y1, term2); x_d2 <= x_d1; v2 <= v1; end
    end

    // STAGE 3
    wire signed [31:0] term3 = TWO - qmult(x_d2, y2);
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin y3 <= 0; x_d3 <= 0; v3 <= 0; end 
        else begin y3 <= qmult(y2, term3); x_d3 <= x_d2; v3 <= v2; end
    end

    // STAGE 4
    wire signed [31:0] term4 = TWO - qmult(x_d3, y3);
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin y4 <= 0; x_d4 <= 0; v4 <= 0; end 
        else begin y4 <= qmult(y3, term4); x_d4 <= x_d3; v4 <= v3; end
    end

    // STAGE 5
    wire signed [31:0] term5 = TWO - qmult(x_d4, y4);
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin y5 <= 0; x_d5 <= 0; v5 <= 0; end 
        else begin y5 <= qmult(y4, term5); x_d5 <= x_d4; v5 <= v4; end
    end

    // STAGE 6 (Final Output) 
    wire signed [31:0] term6 = TWO - qmult(x_d5, y5);
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin y_out <= 0; valid_out <= 0; end
        else if (v5) begin y_out <= qmult(y5, term6); valid_out <= 1; end
        else valid_out <= 0;
    end

endmodule