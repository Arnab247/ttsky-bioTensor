`default_nettype none

module atom_regfile (
    input  wire        clk,
    input  wire        rst_n,
    
    // Write Port
    input  wire        we,
    input  wire [5:0]  w_addr,
    input  wire signed [31:0] w_x, w_y, w_z,
    input  wire [4:0]  w_res_id,   
    input  wire [3:0]  w_atom_idx, 
    
    // Read Ports (4-atom sliding window)
    input  wire [5:0]  r_addr_a, r_addr_b, r_addr_c, r_addr_d,
    
    output wire signed [31:0] xa, ya, za, xb, yb, zb, xc, yc, zc, xd, yd, zd,
    
    // Identity Outputs 
    output wire [4:0]  ra_res_id,  
    output wire [3:0]  ra_atom_idx,
    
    // NEW: Identity outputs for the interacting atom (Atom B)
    output wire [4:0]  rb_res_id,  
    output wire [3:0]  rb_atom_idx 
);

    // Memory arrays for Cartesian Coordinates
    reg signed [31:0] mem_x [0:63];
    reg signed [31:0] mem_y [0:63];
    reg signed [31:0] mem_z [0:63];

    // Memory arrays for Identity Data
    reg [4:0] mem_res_id   [0:63];
    reg [3:0] mem_atom_idx [0:63];

    integer i;

    // Synchronous Write
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            for (i = 0; i < 64; i = i + 1) begin
                mem_x[i] <= 0; mem_y[i] <= 0; mem_z[i] <= 0;
                mem_res_id[i] <= 0; mem_atom_idx[i] <= 0;
            end
        end else if (we) begin
            mem_x[w_addr] <= w_x;
            mem_y[w_addr] <= w_y;
            mem_z[w_addr] <= w_z;
            mem_res_id[w_addr]   <= w_res_id;
            mem_atom_idx[w_addr] <= w_atom_idx;
        end
    end

    // Asynchronous Read (Combinational)
    assign xa = mem_x[r_addr_a]; assign ya = mem_y[r_addr_a]; assign za = mem_z[r_addr_a];
    assign xb = mem_x[r_addr_b]; assign yb = mem_y[r_addr_b]; assign zb = mem_z[r_addr_b];
    assign xc = mem_x[r_addr_c]; assign yc = mem_y[r_addr_c]; assign zc = mem_z[r_addr_c];
    assign xd = mem_x[r_addr_d]; assign yd = mem_y[r_addr_d]; assign zd = mem_z[r_addr_d];

    // Identity Output for Atom A
    assign ra_res_id   = mem_res_id[r_addr_a];
    assign ra_atom_idx = mem_atom_idx[r_addr_a];

    // NEW: Identity Output for Atom B
    assign rb_res_id   = mem_res_id[r_addr_b];
    assign rb_atom_idx = mem_atom_idx[r_addr_b];

endmodule