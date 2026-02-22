`default_nettype none

module residue_database (
    input  wire [4:0]  residue_id,   // e.g., 0=GLY, 1=ALA
    input  wire [3:0]  atom_index,   // Which atom (0..N)
    output reg  [31:0] param_out,    // Mass, Charge, Type
    output reg         is_last_atom, // Flag: End of residue?
    output reg         is_backbone   // Flag: Is this N, CA, or C?
);

    always @(*) begin
        // Default values to prevent latches
        is_last_atom = 0;
        is_backbone  = 0;
        param_out    = 0;

        case (residue_id)
                        // --- 0. GLY ---
            5'd0: begin
                case (atom_index)
                    4'd0: begin param_out = {8'd14, 8'hE2, 16'h0010}; is_backbone = 1; end // N (NH1)
                    4'd1: begin param_out = {8'd12, 8'hFF, 16'h000A}; is_backbone = 1; end // CA (CT2)
                    4'd2: begin param_out = {8'd12, 8'h21, 16'h0001}; is_backbone = 1; end // C (C)
                    4'd3: begin param_out = {8'd16, 8'hDF, 16'h0014}; is_backbone = 1; is_last_atom = 1; end // O (O)
                    default: param_out = 0;
                endcase
            end

            // --- 1. ALA ---
            5'd1: begin
                case (atom_index)
                    4'd0: begin param_out = {8'd14, 8'hE2, 16'h0010}; is_backbone = 1; end // N (NH1)
                    4'd1: begin param_out = {8'd12, 8'h04, 16'h0009}; is_backbone = 1; end // CA (CT1)
                    4'd2: begin param_out = {8'd12, 8'hEF, 16'h000C}; is_backbone = 0; end // CB (CT3)
                    4'd3: begin param_out = {8'd12, 8'h21, 16'h0001}; is_backbone = 1; end // C (C)
                    4'd4: begin param_out = {8'd16, 8'hDF, 16'h0014}; is_backbone = 1; is_last_atom = 1; end // O (O)
                    default: param_out = 0;
                endcase
            end

            // --- 2. VAL ---
            5'd2: begin
                case (atom_index)
                    4'd0: begin param_out = {8'd14, 8'hE2, 16'h0010}; is_backbone = 1; end // N (NH1)
                    4'd1: begin param_out = {8'd12, 8'h04, 16'h0009}; is_backbone = 1; end // CA (CT1)
                    4'd2: begin param_out = {8'd12, 8'hFA, 16'h0009}; is_backbone = 0; end // CB (CT1)
                    4'd3: begin param_out = {8'd12, 8'hEF, 16'h000C}; is_backbone = 0; end // CG1 (CT3)
                    4'd4: begin param_out = {8'd12, 8'hEF, 16'h000C}; is_backbone = 0; end // CG2 (CT3)
                    4'd5: begin param_out = {8'd12, 8'h21, 16'h0001}; is_backbone = 1; end // C (C)
                    4'd6: begin param_out = {8'd16, 8'hDF, 16'h0014}; is_backbone = 1; is_last_atom = 1; end // O (O)
                    default: param_out = 0;
                endcase
            end

            // --- 3. LEU ---
            5'd3: begin
                case (atom_index)
                    4'd0: begin param_out = {8'd14, 8'hE2, 16'h0010}; is_backbone = 1; end // N (NH1)
                    4'd1: begin param_out = {8'd12, 8'h04, 16'h0009}; is_backbone = 1; end // CA (CT1)
                    4'd2: begin param_out = {8'd12, 8'hF4, 16'h000A}; is_backbone = 0; end // CB (CT2)
                    4'd3: begin param_out = {8'd12, 8'hFA, 16'h0009}; is_backbone = 0; end // CG (CT1)
                    4'd4: begin param_out = {8'd12, 8'hEF, 16'h000C}; is_backbone = 0; end // CD1 (CT3)
                    4'd5: begin param_out = {8'd12, 8'hEF, 16'h000C}; is_backbone = 0; end // CD2 (CT3)
                    4'd6: begin param_out = {8'd12, 8'h21, 16'h0001}; is_backbone = 1; end // C (C)
                    4'd7: begin param_out = {8'd16, 8'hDF, 16'h0014}; is_backbone = 1; is_last_atom = 1; end // O (O)
                    default: param_out = 0;
                endcase
            end

            // --- 4. ILE ---
            5'd4: begin
                case (atom_index)
                    4'd0: begin param_out = {8'd14, 8'hE2, 16'h0010}; is_backbone = 1; end // N (NH1)
                    4'd1: begin param_out = {8'd12, 8'h04, 16'h0009}; is_backbone = 1; end // CA (CT1)
                    4'd2: begin param_out = {8'd12, 8'hFA, 16'h0009}; is_backbone = 0; end // CB (CT1)
                    4'd3: begin param_out = {8'd12, 8'hEF, 16'h000C}; is_backbone = 0; end // CG2 (CT3)
                    4'd4: begin param_out = {8'd12, 8'hF4, 16'h000A}; is_backbone = 0; end // CG1 (CT2)
                    4'd5: begin param_out = {8'd12, 8'hEF, 16'h000C}; is_backbone = 0; end // CD (CT3)
                    4'd6: begin param_out = {8'd12, 8'h21, 16'h0001}; is_backbone = 1; end // C (C)
                    4'd7: begin param_out = {8'd16, 8'hDF, 16'h0014}; is_backbone = 1; is_last_atom = 1; end // O (O)
                    default: param_out = 0;
                endcase
            end

            // --- 5. SER ---
            5'd5: begin
                case (atom_index)
                    4'd0: begin param_out = {8'd14, 8'hE2, 16'h0010}; is_backbone = 1; end // N (NH1)
                    4'd1: begin param_out = {8'd12, 8'h04, 16'h0009}; is_backbone = 1; end // CA (CT1)
                    4'd2: begin param_out = {8'd12, 8'h03, 16'h000A}; is_backbone = 0; end // CB (CT2)
                    4'd3: begin param_out = {8'd16, 8'hD6, 16'h0016}; is_backbone = 0; end // OG (OH1)
                    4'd4: begin param_out = {8'd12, 8'h21, 16'h0001}; is_backbone = 1; end // C (C)
                    4'd5: begin param_out = {8'd16, 8'hDF, 16'h0014}; is_backbone = 1; is_last_atom = 1; end // O (O)
                    default: param_out = 0;
                endcase
            end

            // --- 6. THR ---
            5'd6: begin
                case (atom_index)
                    4'd0: begin param_out = {8'd14, 8'hE2, 16'h0010}; is_backbone = 1; end // N (NH1)
                    4'd1: begin param_out = {8'd12, 8'h04, 16'h0009}; is_backbone = 1; end // CA (CT1)
                    4'd2: begin param_out = {8'd12, 8'h09, 16'h0009}; is_backbone = 0; end // CB (CT1)
                    4'd3: begin param_out = {8'd16, 8'hD6, 16'h0016}; is_backbone = 0; end // OG1 (OH1)
                    4'd4: begin param_out = {8'd12, 8'hEF, 16'h000C}; is_backbone = 0; end // CG2 (CT3)
                    4'd5: begin param_out = {8'd12, 8'h21, 16'h0001}; is_backbone = 1; end // C (C)
                    4'd6: begin param_out = {8'd16, 8'hDF, 16'h0014}; is_backbone = 1; is_last_atom = 1; end // O (O)
                    default: param_out = 0;
                endcase
            end

            // --- 7. CYS ---
            5'd7: begin
                case (atom_index)
                    4'd0: begin param_out = {8'd14, 8'hE2, 16'h0010}; is_backbone = 1; end // N (NH1)
                    4'd1: begin param_out = {8'd12, 8'h04, 16'h0009}; is_backbone = 1; end // CA (CT1)
                    4'd2: begin param_out = {8'd12, 8'hF9, 16'h000A}; is_backbone = 0; end // CB (CT2)
                    4'd3: begin param_out = {8'd32, 8'hF1, 16'h0017}; is_backbone = 0; end // SG (S)
                    4'd4: begin param_out = {8'd12, 8'h21, 16'h0001}; is_backbone = 1; end // C (C)
                    4'd5: begin param_out = {8'd16, 8'hDF, 16'h0014}; is_backbone = 1; is_last_atom = 1; end // O (O)
                    default: param_out = 0;
                endcase
            end

            // --- 8. MET ---
            5'd8: begin
                case (atom_index)
                    4'd0: begin param_out = {8'd14, 8'hE2, 16'h0010}; is_backbone = 1; end // N (NH1)
                    4'd1: begin param_out = {8'd12, 8'h04, 16'h0009}; is_backbone = 1; end // CA (CT1)
                    4'd2: begin param_out = {8'd12, 8'hF4, 16'h000A}; is_backbone = 0; end // CB (CT2)
                    4'd3: begin param_out = {8'd12, 8'hF7, 16'h000A}; is_backbone = 0; end // CG (CT2)
                    4'd4: begin param_out = {8'd32, 8'hFA, 16'h0017}; is_backbone = 0; end // SD (S)
                    4'd5: begin param_out = {8'd12, 8'hF2, 16'h000C}; is_backbone = 0; end // CE (CT3)
                    4'd6: begin param_out = {8'd12, 8'h21, 16'h0001}; is_backbone = 1; end // C (C)
                    4'd7: begin param_out = {8'd16, 8'hDF, 16'h0014}; is_backbone = 1; is_last_atom = 1; end // O (O)
                    default: param_out = 0;
                endcase
            end

            // --- 9. ASP ---
            5'd9: begin
                case (atom_index)
                    4'd0: begin param_out = {8'd14, 8'hE2, 16'h0010}; is_backbone = 1; end // N (NH1)
                    4'd1: begin param_out = {8'd12, 8'h04, 16'h0009}; is_backbone = 1; end // CA (CT1)
                    4'd2: begin param_out = {8'd12, 8'hEE, 16'h000B}; is_backbone = 0; end // CB (CT2A)
                    4'd3: begin param_out = {8'd12, 8'h28, 16'h0004}; is_backbone = 0; end // CG (CC)
                    4'd4: begin param_out = {8'd16, 8'hCF, 16'h0015}; is_backbone = 0; end // OD1 (OC)
                    4'd5: begin param_out = {8'd16, 8'hCF, 16'h0015}; is_backbone = 0; end // OD2 (OC)
                    4'd6: begin param_out = {8'd12, 8'h21, 16'h0001}; is_backbone = 1; end // C (C)
                    4'd7: begin param_out = {8'd16, 8'hDF, 16'h0014}; is_backbone = 1; is_last_atom = 1; end // O (O)
                    default: param_out = 0;
                endcase
            end

            // --- 10. GLU ---
            5'd10: begin
                case (atom_index)
                    4'd0: begin param_out = {8'd14, 8'hE2, 16'h0010}; is_backbone = 1; end // N (NH1)
                    4'd1: begin param_out = {8'd12, 8'h04, 16'h0009}; is_backbone = 1; end // CA (CT1)
                    4'd2: begin param_out = {8'd12, 8'hF4, 16'h000B}; is_backbone = 0; end // CB (CT2A)
                    4'd3: begin param_out = {8'd12, 8'hEE, 16'h000A}; is_backbone = 0; end // CG (CT2)
                    4'd4: begin param_out = {8'd12, 8'h28, 16'h0004}; is_backbone = 0; end // CD (CC)
                    4'd5: begin param_out = {8'd16, 8'hCF, 16'h0015}; is_backbone = 0; end // OE1 (OC)
                    4'd6: begin param_out = {8'd16, 8'hCF, 16'h0015}; is_backbone = 0; end // OE2 (OC)
                    4'd7: begin param_out = {8'd12, 8'h21, 16'h0001}; is_backbone = 1; end // C (C)
                    4'd8: begin param_out = {8'd16, 8'hDF, 16'h0014}; is_backbone = 1; is_last_atom = 1; end // O (O)
                    default: param_out = 0;
                endcase
            end

            // --- 11. ASN ---
            5'd11: begin
                case (atom_index)
                    4'd0: begin param_out = {8'd14, 8'hE2, 16'h0010}; is_backbone = 1; end // N (NH1)
                    4'd1: begin param_out = {8'd12, 8'h04, 16'h0009}; is_backbone = 1; end // CA (CT1)
                    4'd2: begin param_out = {8'd12, 8'hF4, 16'h000A}; is_backbone = 0; end // CB (CT2)
                    4'd3: begin param_out = {8'd12, 8'h23, 16'h0004}; is_backbone = 0; end // CG (CC)
                    4'd4: begin param_out = {8'd16, 8'hDD, 16'h0014}; is_backbone = 0; end // OD1 (O)
                    4'd5: begin param_out = {8'd14, 8'hD8, 16'h0011}; is_backbone = 0; end // ND2 (NH2)
                    4'd6: begin param_out = {8'd12, 8'h21, 16'h0001}; is_backbone = 1; end // C (C)
                    4'd7: begin param_out = {8'd16, 8'hDF, 16'h0014}; is_backbone = 1; is_last_atom = 1; end // O (O)
                    default: param_out = 0;
                endcase
            end

            // --- 12. GLN ---
            5'd12: begin
                case (atom_index)
                    4'd0: begin param_out = {8'd14, 8'hE2, 16'h0010}; is_backbone = 1; end // N (NH1)
                    4'd1: begin param_out = {8'd12, 8'h04, 16'h0009}; is_backbone = 1; end // CA (CT1)
                    4'd2: begin param_out = {8'd12, 8'hF4, 16'h000A}; is_backbone = 0; end // CB (CT2)
                    4'd3: begin param_out = {8'd12, 8'hF4, 16'h000A}; is_backbone = 0; end // CG (CT2)
                    4'd4: begin param_out = {8'd12, 8'h23, 16'h0004}; is_backbone = 0; end // CD (CC)
                    4'd5: begin param_out = {8'd16, 8'hDD, 16'h0014}; is_backbone = 0; end // OE1 (O)
                    4'd6: begin param_out = {8'd14, 8'hD8, 16'h0011}; is_backbone = 0; end // NE2 (NH2)
                    4'd7: begin param_out = {8'd12, 8'h21, 16'h0001}; is_backbone = 1; end // C (C)
                    4'd8: begin param_out = {8'd16, 8'hDF, 16'h0014}; is_backbone = 1; is_last_atom = 1; end // O (O)
                    default: param_out = 0;
                endcase
            end

            // --- 13. LYS ---
            5'd13: begin
                case (atom_index)
                    4'd0: begin param_out = {8'd14, 8'hE2, 16'h0010}; is_backbone = 1; end // N (NH1)
                    4'd1: begin param_out = {8'd12, 8'h04, 16'h0009}; is_backbone = 1; end // CA (CT1)
                    4'd2: begin param_out = {8'd12, 8'hF4, 16'h000A}; is_backbone = 0; end // CB (CT2)
                    4'd3: begin param_out = {8'd12, 8'hF4, 16'h000A}; is_backbone = 0; end // CG (CT2)
                    4'd4: begin param_out = {8'd12, 8'hF4, 16'h000A}; is_backbone = 0; end // CD (CT2)
                    4'd5: begin param_out = {8'd12, 8'h0D, 16'h000A}; is_backbone = 0; end // CE (CT2)
                    4'd6: begin param_out = {8'd14, 8'hED, 16'h0012}; is_backbone = 0; end // NZ (NH3)
                    4'd7: begin param_out = {8'd12, 8'h21, 16'h0001}; is_backbone = 1; end // C (C)
                    4'd8: begin param_out = {8'd16, 8'hDF, 16'h0014}; is_backbone = 1; is_last_atom = 1; end // O (O)
                    default: param_out = 0;
                endcase
            end

            // --- 14. ARG ---
            5'd14: begin
                case (atom_index)
                    4'd0: begin param_out = {8'd14, 8'hE2, 16'h0010}; is_backbone = 1; end // N (NH1)
                    4'd1: begin param_out = {8'd12, 8'h04, 16'h0009}; is_backbone = 1; end // CA (CT1)
                    4'd2: begin param_out = {8'd12, 8'hF4, 16'h000A}; is_backbone = 0; end // CB (CT2)
                    4'd3: begin param_out = {8'd12, 8'hF4, 16'h000A}; is_backbone = 0; end // CG (CT2)
                    4'd4: begin param_out = {8'd12, 8'h0D, 16'h000A}; is_backbone = 0; end // CD (CT2)
                    4'd5: begin param_out = {8'd14, 8'hD3, 16'h000F}; is_backbone = 0; end // NE (NC2)
                    4'd6: begin param_out = {8'd12, 8'h29, 16'h0001}; is_backbone = 0; end // CZ (C)
                    4'd7: begin param_out = {8'd14, 8'hCD, 16'h000F}; is_backbone = 0; end // NH1 (NC2)
                    4'd8: begin param_out = {8'd14, 8'hCD, 16'h000F}; is_backbone = 0; end // NH2 (NC2)
                    4'd9: begin param_out = {8'd12, 8'h21, 16'h0001}; is_backbone = 1; end // C (C)
                    4'd10: begin param_out = {8'd16, 8'hDF, 16'h0014}; is_backbone = 1; is_last_atom = 1; end // O (O)
                    default: param_out = 0;
                endcase
            end

            // --- 15. PHE ---
            5'd15: begin
                case (atom_index)
                    4'd0: begin param_out = {8'd14, 8'hE2, 16'h0010}; is_backbone = 1; end // N (NH1)
                    4'd1: begin param_out = {8'd12, 8'h04, 16'h0009}; is_backbone = 1; end // CA (CT1)
                    4'd2: begin param_out = {8'd12, 8'hF4, 16'h000A}; is_backbone = 0; end // CB (CT2)
                    4'd3: begin param_out = {8'd12, 8'h00, 16'h0002}; is_backbone = 0; end // CG (CA)
                    4'd4: begin param_out = {8'd12, 8'hF9, 16'h0002}; is_backbone = 0; end // CD1 (CA)
                    4'd5: begin param_out = {8'd12, 8'hF9, 16'h0002}; is_backbone = 0; end // CE1 (CA)
                    4'd6: begin param_out = {8'd12, 8'hF9, 16'h0002}; is_backbone = 0; end // CZ (CA)
                    4'd7: begin param_out = {8'd12, 8'hF9, 16'h0002}; is_backbone = 0; end // CD2 (CA)
                    4'd8: begin param_out = {8'd12, 8'hF9, 16'h0002}; is_backbone = 0; end // CE2 (CA)
                    4'd9: begin param_out = {8'd12, 8'h21, 16'h0001}; is_backbone = 1; end // C (C)
                    4'd10: begin param_out = {8'd16, 8'hDF, 16'h0014}; is_backbone = 1; is_last_atom = 1; end // O (O)
                    default: param_out = 0;
                endcase
            end

            // --- 16. TYR ---
            5'd16: begin
                case (atom_index)
                    4'd0: begin param_out = {8'd14, 8'hE2, 16'h0010}; is_backbone = 1; end // N (NH1)
                    4'd1: begin param_out = {8'd12, 8'h04, 16'h0009}; is_backbone = 1; end // CA (CT1)
                    4'd2: begin param_out = {8'd12, 8'hF4, 16'h000A}; is_backbone = 0; end // CB (CT2)
                    4'd3: begin param_out = {8'd12, 8'h00, 16'h0002}; is_backbone = 0; end // CG (CA)
                    4'd4: begin param_out = {8'd12, 8'hF9, 16'h0002}; is_backbone = 0; end // CD1 (CA)
                    4'd5: begin param_out = {8'd12, 8'hF9, 16'h0002}; is_backbone = 0; end // CE1 (CA)
                    4'd6: begin param_out = {8'd12, 8'h07, 16'h0002}; is_backbone = 0; end // CZ (CA)
                    4'd7: begin param_out = {8'd16, 8'hDD, 16'h0016}; is_backbone = 0; end // OH (OH1)
                    4'd8: begin param_out = {8'd12, 8'hF9, 16'h0002}; is_backbone = 0; end // CD2 (CA)
                    4'd9: begin param_out = {8'd12, 8'hF9, 16'h0002}; is_backbone = 0; end // CE2 (CA)
                    4'd10: begin param_out = {8'd12, 8'h21, 16'h0001}; is_backbone = 1; end // C (C)
                    4'd11: begin param_out = {8'd16, 8'hDF, 16'h0014}; is_backbone = 1; is_last_atom = 1; end // O (O)
                    default: param_out = 0;
                endcase
            end

            // --- 17: HISTIDINE (HIS/HSD) CORRECTED ---
            5'd17: begin
                case (atom_index)
                    4'd0: begin // N (Charge: -0.47)
                        param_out = {8'd14, 8'hE2, 16'h0014};
                        is_backbone = 1;
                    end
                    4'd1: begin // CA (Charge: 0.07)
                        param_out = {8'd12, 8'h04, 16'h0015};
                        is_backbone = 1;
                    end
                    4'd2: begin // CB (Charge: -0.09)
                        param_out = {8'd12, 8'hFA, 16'h0016};
                        is_backbone = 0;
                    end
                    4'd3: begin // ND1 (Charge: -0.36)
                        param_out = {8'd14, 8'hE9, 16'h0017};
                        is_backbone = 0;
                    end
                    4'd4: begin // CG (Charge: -0.05)
                        param_out = {8'd12, 8'hFD, 16'h0018};
                        is_backbone = 0;
                    end
                    4'd5: begin // CE1 (Charge: 0.25)
                        param_out = {8'd12, 8'h10, 16'h0019};
                        is_backbone = 0;
                    end
                    4'd6: begin // NE2 (Charge: -0.7)
                        param_out = {8'd14, 8'hD3, 16'h001A};
                        is_backbone = 0;
                    end
                    4'd7: begin // CD2 (Charge: 0.22)
                        param_out = {8'd12, 8'h0E, 16'h001B};
                        is_backbone = 0;
                    end
                    4'd8: begin // C (Charge: 0.51)
                        param_out = {8'd12, 8'h21, 16'h001C};
                        is_backbone = 1;
                    end
                    4'd9: begin // O (Charge: -0.51)
                        param_out = {8'd16, 8'hDF, 16'h001D};
                        is_backbone = 0;
                        is_last_atom = 1;
                    end
                    default: begin
                        param_out = 0;
                        is_last_atom = 1;
                    end
                endcase
            end

            // --- 18. TRP ---
            5'd18: begin
                case (atom_index)
                    4'd0: begin param_out = {8'd14, 8'hE2, 16'h0010}; is_backbone = 1; end // N (NH1)
                    4'd1: begin param_out = {8'd12, 8'h04, 16'h0009}; is_backbone = 1; end // CA (CT1)
                    4'd2: begin param_out = {8'd12, 8'hF4, 16'h000A}; is_backbone = 0; end // CB (CT2)
                    4'd3: begin param_out = {8'd12, 8'hFE, 16'h000D}; is_backbone = 0; end // CG (CY)
                    4'd4: begin param_out = {8'd12, 8'hF6, 16'h0002}; is_backbone = 0; end // CD1 (CA)
                    4'd5: begin param_out = {8'd14, 8'hDF, 16'h0013}; is_backbone = 0; end // NE1 (NY)
                    4'd6: begin param_out = {8'd12, 8'h0F, 16'h0008}; is_backbone = 0; end // CE2 (CPT)
                    4'd7: begin param_out = {8'd12, 8'h07, 16'h0008}; is_backbone = 0; end // CD2 (CPT)
                    4'd8: begin param_out = {8'd12, 8'hF0, 16'h0003}; is_backbone = 0; end // CE3 (CAI)
                    4'd9: begin param_out = {8'd12, 8'hF3, 16'h0002}; is_backbone = 0; end // CZ3 (CA)
                    4'd10: begin param_out = {8'd12, 8'hEF, 16'h0003}; is_backbone = 0; end // CZ2 (CAI)
                    4'd11: begin param_out = {8'd12, 8'hF7, 16'h0002}; is_backbone = 0; end // CH2 (CA)
                    4'd12: begin param_out = {8'd12, 8'h21, 16'h0001}; is_backbone = 1; end // C (C)
                    4'd13: begin param_out = {8'd16, 8'hDF, 16'h0014}; is_backbone = 1; is_last_atom = 1; end // O (O)
                    default: param_out = 0;
                endcase
            end

            // --- 19. PRO ---
            5'd19: begin
                case (atom_index)
                    4'd0: begin param_out = {8'd14, 8'hED, 16'h000E}; is_backbone = 1; end // N (N)
                    4'd1: begin param_out = {8'd12, 8'h00, 16'h0007}; is_backbone = 0; end // CD (CP3)
                    4'd2: begin param_out = {8'd12, 8'h01, 16'h0005}; is_backbone = 1; end // CA (CP1)
                    4'd3: begin param_out = {8'd12, 8'hF4, 16'h0006}; is_backbone = 0; end // CB (CP2)
                    4'd4: begin param_out = {8'd12, 8'hF4, 16'h0006}; is_backbone = 0; end // CG (CP2)
                    4'd5: begin param_out = {8'd12, 8'h21, 16'h0001}; is_backbone = 1; end // C (C)
                    4'd6: begin param_out = {8'd16, 8'hDF, 16'h0014}; is_backbone = 1; is_last_atom = 1; end // O (O)
                    default: param_out = 0;
                endcase
            end

            // --- DEFAULT: Unknown Residue ID ---
            // Note: Residue ID 17 (HIS - Histidine) is missing and reserved for future implementation
            // IDs 20-31 are also unused
            default: begin
                // Return zeros for unknown residues
                param_out = 32'h00000000;
                is_last_atom = 0;
                is_backbone = 0;
            end

        endcase
    end
endmodule