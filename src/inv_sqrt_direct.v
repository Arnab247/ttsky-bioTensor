`default_nettype none

module inv_sqrt_direct (
    input  wire clk,
    input  wire rst_n,
    input  wire valid_in,
    input  wire signed [31:0] x,
    output reg  signed [31:0] y_out,
    output reg  valid_out
);

    // Q16.16 Multiplier Function
    function automatic signed [31:0] qmult(input signed [31:0] a, input signed [31:0] b);
        reg signed [63:0] temp;
        begin
            temp = $signed({{32{a[31]}}, a}) * $signed({{32{b[31]}}, b});
            qmult = temp[47:16];
        end
    endfunction

    // 1.5 in Q16.16
    wire signed [31:0] CONST_1_5 = 32'h0001_8000; 

    // --- PIPELINE REGISTERS ---
    reg [5:0] v_sr;
    reg signed [4:0] k_d1, k_d2, k_d3, k_d4, k_d5; // Shift amount shadow registers
    reg signed [31:0] x_norm_d1, x_norm_d2, x_norm_d3;
    reg signed [31:0] y0_d1, y0_d2, y0_d3;
    
    reg signed [31:0] y0_sq;
    reg signed [31:0] x_y0_sq;
    reg signed [31:0] sub_term;
    reg signed [31:0] y1;

    // --- COMBINATIONAL: Leading Zero Detector (LZD) ---
    reg [3:0] P;
    always @(*) begin
        if      (x[31:30] != 0) P = 15;
        else if (x[29:28] != 0) P = 14;
        else if (x[27:26] != 0) P = 13;
        else if (x[25:24] != 0) P = 12;
        else if (x[23:22] != 0) P = 11;
        else if (x[21:20] != 0) P = 10;
        else if (x[19:18] != 0) P = 9;
        else if (x[17:16] != 0) P = 8;
        else if (x[15:14] != 0) P = 7;
        else if (x[13:12] != 0) P = 6;
        else if (x[11:10] != 0) P = 5;
        else if (x[9:8]   != 0) P = 4;
        else if (x[7:6]   != 0) P = 3;
        else if (x[5:4]   != 0) P = 2;
        else if (x[3:2]   != 0) P = 1;
        else                    P = 0;
    end

    wire signed [4:0] k_shift = 8 - P; 
    wire signed [31:0] x_norm = (k_shift >= 0) ? (x << (2 * k_shift)) : (x >> (-2 * k_shift));

    // --- STAGE 1: Normalize & Seed ---
    reg signed [31:0] y0;
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            y0 <= 0; x_norm_d1 <= 0; k_d1 <= 0;
        end else begin
            x_norm_d1 <= x_norm;
            k_d1 <= k_shift;
            
            // 24-entry LUT for values between 1.0 and 4.0
            case (x_norm[17:13])
                5'd8:  y0 <= 32'h0000_F858; // 0.970
                5'd9:  y0 <= 32'h0000_EAE7; // 0.917
                5'd10: y0 <= 32'h0000_DF7A; // 0.872
                5'd11: y0 <= 32'h0000_D57A; // 0.833
                5'd12: y0 <= 32'h0000_CCCD; // 0.800
                5'd13: y0 <= 32'h0000_C511; // 0.769
                5'd14: y0 <= 32'h0000_BE22; // 0.742
                5'd15: y0 <= 32'h0000_B7E8; // 0.718
                5'd16: y0 <= 32'h0000_B241; // 0.696
                5'd17: y0 <= 32'h0000_AD15; // 0.676
                5'd18: y0 <= 32'h0000_A853; // 0.657
                5'd19: y0 <= 32'h0000_A3F0; // 0.640
                5'd20: y0 <= 32'h0000_9FE9; // 0.624
                5'd21: y0 <= 32'h0000_9C25; // 0.609
                5'd22: y0 <= 32'h0000_98A2; // 0.596
                5'd23: y0 <= 32'h0000_955B; // 0.583
                5'd24: y0 <= 32'h0000_924D; // 0.571
                5'd25: y0 <= 32'h0000_8F6B; // 0.560
                5'd26: y0 <= 32'h0000_8CBA; // 0.549
                5'd27: y0 <= 32'h0000_8A23; // 0.539
                5'd28: y0 <= 32'h0000_87A8; // 0.529
                5'd29: y0 <= 32'h0000_8559; // 0.520
                5'd30: y0 <= 32'h0000_831F; // 0.512
                5'd31: y0 <= 32'h0000_8100; // 0.503
                default: y0 <= 32'h0001_0000;
            endcase
        end
    end

    // --- STAGE 2: Square the guess (y0^2) ---
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin y0_sq <= 0; x_norm_d2 <= 0; y0_d1 <= 0; k_d2 <= 0; end
        else begin
            y0_sq <= qmult(y0, y0);
            x_norm_d2 <= x_norm_d1;
            y0_d1 <= y0;
            k_d2 <= k_d1;
        end
    end

    // --- STAGE 3: Multiply by x (x_norm * y0^2) ---
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin x_y0_sq <= 0; y0_d2 <= 0; k_d3 <= 0; end
        else begin
            x_y0_sq <= qmult(x_norm_d2, y0_sq);
            y0_d2 <= y0_d1;
            k_d3 <= k_d2;
        end
    end

    // --- STAGE 4: Multiply by 0.5 and Subtract (1.5 - 0.5 * x * y0^2) ---
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin sub_term <= 0; y0_d3 <= 0; k_d4 <= 0; end
        else begin
            sub_term <= CONST_1_5 - (x_y0_sq >>> 1); // >>> 1 is multiplying by 0.5
            y0_d3 <= y0_d2;
            k_d4 <= k_d3;
        end
    end

    // --- STAGE 5: Final Multiply (y0 * (1.5 - ...)) ---
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin y1 <= 0; k_d5 <= 0; end
        else begin
            y1 <= qmult(y0_d3, sub_term);
            k_d5 <= k_d4;
        end
    end

    // --- STAGE 6: Denormalize (Output Result) ---
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) y_out <= 0;
        else begin
            // Shift back by half the amount we shifted initially
            if (k_d5 >= 0) y_out <= (y1 << k_d5);
            else           y_out <= (y1 >> -k_d5);
        end
    end

    // --- VALID SIGNAL PIPELINE ---
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            v_sr <= 6'b0;
            valid_out <= 0;
        end else begin
            v_sr <= {v_sr[4:0], valid_in};
            valid_out <= v_sr[5]; // Aligns with Stage 6
        end
    end

endmodule