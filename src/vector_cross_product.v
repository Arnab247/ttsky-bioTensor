`default_nettype none

module vector_cross_product (
    input  wire signed [31:0] Ax, Ay, Az,
    input  wire signed [31:0] Bx, By, Bz,
    output wire signed [31:0] Rx, Ry, Rz
);

    function automatic signed [31:0] qmult(input signed [31:0] a, input signed [31:0] b);
        reg signed [63:0] temp;
        begin
            temp = $signed({{32{a[31]}}, a}) * $signed({{32{b[31]}}, b});
            qmult = temp[47:16];
        end
    endfunction

    assign Rx = qmult(Ay, Bz) - qmult(Az, By);
    assign Ry = qmult(Az, Bx) - qmult(Ax, Bz);
    assign Rz = qmult(Ax, By) - qmult(Ay, Bx);

endmodule