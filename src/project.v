`default_nettype none

module tt_um_bioTensor (
    input  wire [7:0] ui_in,    // [7:6]=Cmd, [5:0]=Data/Addr
    output wire [7:0] uo_out,   // [7]=Done, [6:0]=Iter Progress
    input  wire [7:0] uio_in,   // Unused
    output wire [7:0] uio_out,  // Unused (must be 0)
    output wire [7:0] uio_oe,   // Unused (must be 0)
    input  wire       ena,      // Power enable [cite: 150]
    input  wire       clk,      // Clock [cite: 149]
    input  wire       rst_n     // Active-low reset [cite: 151]
);

    // --- Mandatory Tiny Tapeout Wiring ---
    assign uio_out = 8'b00000000; // [cite: 156]
    assign uio_oe  = 8'b00000000; // All bidirectional as inputs [cite: 157]
    wire _unused = &{ena, uio_in, 1'b0}; // Prevent warnings [cite: 160]

    // --- Data Sequencer Registers ---
    reg [31:0] shift_reg;
    reg [5:0]  latched_addr;
    reg        load_pulse;

    // --- MD System Wires ---
    wire md_done;
    wire [15:0] current_iter;

    // Command Decoding:
    // 00: Load Byte 0 | 01: Load Byte 1 | 10: Load Byte 2 | 11: Load Byte 3 + Trigger Write
    always @(posedge clk) begin
        if (!rst_n) begin
            shift_reg <= 0;
            load_pulse <= 0;
        end else begin
            load_pulse <= 0;
            case (ui_in[7:6])
                2'b00: shift_reg[7:0]   <= {2'b0, ui_in[5:0]}; // Small coordinates for test
                2'b01: shift_reg[15:8]  <= {2'b0, ui_in[5:0]};
                2'b10: shift_reg[23:16] <= {2'b0, ui_in[5:0]};
                2'b11: begin
                    shift_reg[31:24] <= {2'b0, ui_in[5:0]};
                    load_pulse <= 1; // Pulse high for one cycle to write to regfile
                end
            endcase
        end
    end

    // --- Physical Pin Mapping ---
    assign uo_out[7]   = md_done;         // Pin 7: Success flag [cite: 153]
    assign uo_out[6:0] = current_iter[6:0]; // Pins 0-6: Progress [cite: 153]

    // --- MD Core Instantiation ---
    md_system_top user_project (
        .clk        (clk),
        .rst_n      (rst_n),
        .start_run  (ena), // Starts when chip is enabled [cite: 150]
        .max_iters  (16'd100),
        .num_atoms  (6'd10),
        
        // Sequencer Interface
        .load_en    (load_pulse),
        .load_addr  (ui_in[5:0]), 
        .load_x     (shift_reg), 
        .load_y     (32'h0), // For 3D, you'd need more commands
        .load_z     (32'h0),
        
        .done       (md_done),
        .iter_count (current_iter)
    );

endmodule