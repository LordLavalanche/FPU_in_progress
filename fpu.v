`timescale 1ns / 1ps
module cmp(
    input clk,
    input rst,
    input [31:0] a,
    input [31:0] b,
    output reg [2:0] con
);
    always @(posedge clk or posedge rst) begin
        if (rst) begin
            con <= 3'b000;
        end else begin
            con[2] <= a[31]; // sign bit of a
            con[1] <= b[31]; // sign bit of b
            con[0] <= (a[30:0] == b[30:0]) ? 1'b1 : 1'b0; // exponent+mantissa comparison
        end
    end
endmodule
module addsub(input clk, input rst, input [31:0] a, input [31:0] b, input x, output reg [31:0] out);
    reg [2:0] pipe1;
    reg [22:0] c, d; 
    reg [7:0] e;
    reg sign;
    
    always @(posedge clk or posedge rst) begin
        if (rst) begin
            pipe1 <= 3'b0;
        end else begin
            pipe1 <= {a[31], b[31], a[30:0] > b[30:0]};
        end
    end
    
    always @(posedge clk or posedge rst) begin
        if (rst) begin
            c <= 23'b0;
            d <= 23'b0;
            // Set exponent and sign bits in out to 0 on reset
            out[30:23] <= 8'b0;
            out[31] <= 1'b0;
        end else begin
            case (pipe1)
                3'b000, 3'b110: begin // b >= a
                    out[30:23] <= b[30:23];
                    c <= b[22:0];
                    d <= a[22:0] >> (b[30:23] - a[30:23]);
                    out[31] <= b[31];
                end
                3'b001, 3'b111: begin // a > b
                    out[30:23] <= a[30:23];
                    c <= a[22:0];
                    d <= b[22:0] >> (a[30:23] - b[30:23]);
                    out[31] <= a[31];
                end
                3'b010: begin // a[31]=0, b[31]=1
                    out[30:23] <= a[30:23];
                    c <= a[22:0];
                    d <= b[22:0] >> (a[30:23] - b[30:23]);
                    out[31] <= a[31];
                end
                3'b101: begin // a[31]=1, b[31]=0
                    out[30:23] <= b[30:23];
                    c <= b[22:0];
                    d <= a[22:0] >> (b[30:23] - a[30:23]);
                    out[31] <= b[31];
                end
                default: begin
                    c <= 23'b0;
                    d <= 23'b0;
                    out[30:23] <= 8'b0;
                    out[31] <= 1'b0;
                end
            endcase
        end
    end

    // Adder control signals
    reg [22:0] adder_a, adder_b;
    reg adder_cin;
    wire [22:0] adder_sum;
    wire adder_cout;
    
    // Instantiate Kogge-Stone adder outside always block
    kogge_stone ks_inst(.x(adder_a), .y(adder_b), .sum(adder_sum), .cin(adder_cin), .cout(adder_cout));
    
    // Control logic for adder inputs
    always @(posedge clk or posedge rst) begin
        if (rst) begin
            adder_a <= 23'b0;
            adder_b <= 23'b0;
            adder_cin <= 1'b0;
        end else begin
            case (pipe1)
                3'b000: begin // a[31]=0, b[31]=0, b>=a
                    if (x == 1'b0) begin // addition
                        adder_a <= c;
                        adder_b <= d;
                        adder_cin <= 1'b0;
                    end else begin // subtraction: d - c
                        adder_a <= d;
                        adder_b <= ~c;
                        adder_cin <= 1'b1;
                    end
                end
                3'b001: begin // a[31]=0, b[31]=0, a>b
                    if (x == 1'b0) begin // addition
                        adder_a <= c;
                        adder_b <= d;
                        adder_cin <= 1'b0;
                    end else begin // subtraction: c - d
                        adder_a <= c;
                        adder_b <= ~d;
                        adder_cin <= 1'b1;
                    end
                end
                3'b010: begin // a[31]=0, b[31]=1
                    if (x == 1'b0) begin // addition
                        adder_a <= c;
                        adder_b <= d;
                        adder_cin <= 1'b0;
                    end else begin // subtraction: c - d
                        adder_a <= c;
                        adder_b <= ~d;
                        adder_cin <= 1'b1;
                    end
                end
                3'b100: begin // a[31]=1, b[31]=0
                    if (x == 1'b0) begin // addition
                        adder_a <= c;
                        adder_b <= d;
                        adder_cin <= 1'b0;
                    end else begin // subtraction: d - c
                        adder_a <= d;
                        adder_b <= ~c;
                        adder_cin <= 1'b1;
                    end
                end
                3'b110: begin // a[31]=1, b[31]=1, b>=a
                    if (x == 1'b0) begin // addition
                        adder_a <= c;
                        adder_b <= d;
                        adder_cin <= 1'b0;
                    end else begin // subtraction: d - c
                        adder_a <= d;
                        adder_b <= ~c;
                        adder_cin <= 1'b1;
                    end
                end
                3'b111: begin // a[31]=1, b[31]=1, a>b
                    if (x == 1'b0) begin // addition
                        adder_a <= c;
                        adder_b <= d;
                        adder_cin <= 1'b0;
                    end else begin // subtraction: c - d
                        adder_a <= c;
                        adder_b <= ~d;
                        adder_cin <= 1'b1;
                    end
                end
                default: begin
                    adder_a <= 23'b0;
                    adder_b <= 23'b0;
                    adder_cin <= 1'b0;
                end
            endcase
        end
    end
    
    // Assign adder output to result
    always @(posedge clk or posedge rst) begin
        if (rst) begin
            out[22:0] <= 23'b0;
        end else begin
            out[22:0] <= adder_sum;
        end
    end
endmodule


module gray_cell(Gkj, Pik, Gik, G);
    input Gkj, Pik, Gik;
    output G;
    wire Y;
    and(Y, Gkj, Pik);
    or(G, Y, Gik);
endmodule
module black_cell(Gkj, Pik, Gik, Pkj, G, P);
    input Gkj, Pik, Gik, Pkj;
    output G, P;
    wire Y;
    and(Y, Gkj, Pik);
    or(G, Gik, Y);
    and(P, Pkj, Pik);
endmodule
module and_xor(a, b, p, g);
    input a, b;
    output p, g;
    xor(p, a, b);
    and(g, a, b);
endmodule
module kogge_stone(
    input        [22:0]  x,
    input        [22:0]  y,
    input                cin,
    output       [22:0]  sum,
    output               cout
);

    // Number of stages for the prefix network
    localparam STAGES = 5; // ceil(log2(23))

    // Wires for Propagate and Generate signals at each level
    wire [22:0] p [0:STAGES];
    wire [22:0] g [0:STAGES];
    
    // Wire for the final carry into each bit position
    wire [23:0] c;

    // --- Level 0: Initial P and G Calculation ---
    // This stage computes the propagate and generate for each bit.
    genvar i;
    generate
        for (i = 0; i < 23; i = i + 1) begin : precomputation
            and_xor init_pg (
                .a(x[i]), 
                .b(y[i]), 
                .p(p[0][i]), 
                .g(g[0][i])
            );
        end
    endgenerate


    // --- Levels 1 to 5: Parallel Prefix Network ---
    // This is the core of the Kogge-Stone adder. Each level doubles
    // the range of the group P and G signals.
    genvar s; // s for stage
    generate
        for (s = 0; s < STAGES-1; s = s + 1) begin : prefix_network_stage
            genvar j; // j for bit index
            for (j = 0; j < 23; j = j + 1) begin : bit_logic
                
                // If the bit is within the range to have a predecessor
                // at distance 2^s, instantiate a black cell.
                if (j >= (1 << s)) begin
                    black_cell bc (
                        .Gkj(g[s][j - (1 << s)]), // Generate from predecessor
                        .Pkj(p[s][j - (1 << s)]), // Propagate from predecessor
                        .Gik(g[s][j]),           // Generate from current bit
                        .Pik(p[s][j]),           // Propagate from current bit
                        .G(g[s+1][j]),           // Output new Group Generate
                        .P(p[s+1][j])            // Output new Group Propagate
                    );
                end else begin
                    // Otherwise, just pass the signals to the next level.
                    assign g[s+1][j] = g[s][j];
                    assign p[s+1][j] = p[s][j];
                end
                
            end
        end
    endgenerate


    // --- Final Sum and Carry Out Calculation ---

    // The carry-in to bit 0 is the module's carry-in.
    assign c[0] = cin;

    // The carry-in to bit 'i' is the final group generate signal from bit 'i-1'.
    assign c[23:1] = g[STAGES-1][22:0];
    
    // The final carry-out of the adder.
    assign cout = c[23];

    // The sum is the initial propagate XORed with the carry-in for each bit.
    assign sum = p[0] ^ c[22:0];

endmodule
`timescale 1ns / 1ps

// Corrected Floating Point Multiplication Module
module muldiv(
    input               clk,
    input               rst,
    input       [31:0]  a,
    input       [31:0]  b,
    output reg  [31:0]  result
);

    // --- STAGE 0: Combinational Logic (Unpacking & Exponent Math) ---

    // Unpack inputs using continuous assign statements for clarity
    wire sign_a = a[31];
    wire sign_b = b[31];
    wire [7:0] exp_a = a[30:23];
    wire [7:0] exp_b = b[30:23];
    wire [23:0] mant_a = {1'b1, a[22:0]}; // Add implicit 1
    wire [23:0] mant_b = {1'b1, b[22:0]};

    // Combinational logic for exponent calculation
    wire [7:0] exp_sum_ks;
    wire [7:0] exp_result_logic;
    wire [7:0] bias = 8'd127;
    wire [7:0] bias_inv = ~bias;

    // Instantiate adders structurally (outside of always blocks)
    kogge_stone8 exp_adder (
        .x(exp_a),
        .y(exp_b),
        .sum(exp_sum_ks),
        .cin(1'b0),
        .cout() 
    );

    kogge_stone8 bias_subtractor (
        .x(exp_sum_ks),
        .y(bias_inv),
        .sum(exp_result_logic),
        .cin(1'b1), // Add 1 for two's complement subtraction
        .cout()
    );

    // --- Pipeline Registers ---
    // Stage 1 registers (outputs of the first pipeline stage)
    reg         p1_sign;
    reg [7:0]   p1_exp;
    reg [23:0]  p1_mant_a;
    reg [23:0]  p1_mant_b;

    // --- PIPELINE STAGE 1: Latch Unpacked Data ---
    // This stage captures all necessary data on the clock edge. This ensures
    // the mantissas are synchronized with their corresponding sign and exponent.
    always @(posedge clk or posedge rst) begin
        if (rst) begin
            p1_sign   <= 1'b0;
            p1_exp    <= 8'd0;
            p1_mant_a <= 24'd0;
            p1_mant_b <= 24'd0;
        end else begin
            p1_sign   <= sign_a ^ sign_b;
            p1_exp    <= exp_result_logic;
            p1_mant_a <= mant_a; // Latch the mantissa for the next stage
            p1_mant_b <= mant_b; // Latch the mantissa for the next stage
        end
    end

    // --- STAGE 2: Combinational Mantissa Multiplication ---
    // The multiplier now correctly operates on the pipelined data from Stage 1.
    wire [47:0] mant_prod;
    mult wallace_mult (
        .a(p1_mant_a),
        .b(p1_mant_b),
        .product(mant_prod)
    );

    // Wire for incremented exponent
    wire [7:0] exp_incremented = p1_exp + 8'd1;
    wire [7:0] exp_decremented = p1_exp - 8'd1;
    
    // Leading one detection and normalization logic
    reg [4:0] leading_zeros;
    reg [22:0] normalized_mantissa;
    reg [7:0] final_exponent;
    
    // Leading one detection (priority encoder)
    always @(*) begin
        // Default values
        leading_zeros = 5'd0;
        normalized_mantissa = 23'b0;
        final_exponent = p1_exp;
        
        casez (mant_prod[47:23]) // Check upper 25 bits for leading 1
            25'b1????????????????????????: begin // Leading 1 at bit 47
                leading_zeros = 5'd0;
                normalized_mantissa = mant_prod[46:24];
                final_exponent = exp_incremented;
            end
            25'b01???????????????????????: begin // Leading 1 at bit 46
                leading_zeros = 5'd0;
                normalized_mantissa = mant_prod[45:23];
                final_exponent = p1_exp;
            end
            25'b001??????????????????????: begin // Leading 1 at bit 45
                leading_zeros = 5'd1;
                normalized_mantissa = mant_prod[44:22];
                final_exponent = exp_decremented;
            end
            25'b0001?????????????????????: begin // Leading 1 at bit 44
                leading_zeros = 5'd2;
                normalized_mantissa = mant_prod[43:21];
                final_exponent = p1_exp - 8'd2;
            end
            25'b00001????????????????????: begin // Leading 1 at bit 43
                leading_zeros = 5'd3;
                normalized_mantissa = mant_prod[42:20];
                final_exponent = p1_exp - 8'd3;
            end
            25'b000001???????????????????: begin // Leading 1 at bit 42
                leading_zeros = 5'd4;
                normalized_mantissa = mant_prod[41:19];
                final_exponent = p1_exp - 8'd4;
            end
            25'b0000001??????????????????: begin // Leading 1 at bit 41
                leading_zeros = 5'd5;
                normalized_mantissa = mant_prod[40:18];
                final_exponent = p1_exp - 8'd5;
            end
            25'b00000001?????????????????: begin // Leading 1 at bit 40
                leading_zeros = 5'd6;
                normalized_mantissa = mant_prod[39:17];
                final_exponent = p1_exp - 8'd6;
            end
            25'b000000001????????????????: begin // Leading 1 at bit 39
                leading_zeros = 5'd7;
                normalized_mantissa = mant_prod[38:16];
                final_exponent = p1_exp - 8'd7;
            end
            25'b0000000001???????????????: begin // Leading 1 at bit 38
                leading_zeros = 5'd8;
                normalized_mantissa = mant_prod[37:15];
                final_exponent = p1_exp - 8'd8;
            end
            25'b00000000001??????????????: begin // Leading 1 at bit 37
                leading_zeros = 5'd9;
                normalized_mantissa = mant_prod[36:14];
                final_exponent = p1_exp - 8'd9;
            end
            25'b000000000001?????????????: begin // Leading 1 at bit 36
                leading_zeros = 5'd10;
                normalized_mantissa = mant_prod[35:13];
                final_exponent = p1_exp - 8'd10;
            end
            25'b0000000000001????????????: begin // Leading 1 at bit 35
                leading_zeros = 5'd11;
                normalized_mantissa = mant_prod[34:12];
                final_exponent = p1_exp - 8'd11;
            end
            25'b00000000000001???????????: begin // Leading 1 at bit 34
                leading_zeros = 5'd12;
                normalized_mantissa = mant_prod[33:11];
                final_exponent = p1_exp - 8'd12;
            end
            25'b000000000000001??????????: begin // Leading 1 at bit 33
                leading_zeros = 5'd13;
                normalized_mantissa = mant_prod[32:10];
                final_exponent = p1_exp - 8'd13;
            end
            25'b0000000000000001?????????: begin // Leading 1 at bit 32
                leading_zeros = 5'd14;
                normalized_mantissa = mant_prod[31:9];
                final_exponent = p1_exp - 8'd14;
            end
            25'b00000000000000001????????: begin // Leading 1 at bit 31
                leading_zeros = 5'd15;
                normalized_mantissa = mant_prod[30:8];
                final_exponent = p1_exp - 8'd15;
            end
            25'b000000000000000001???????: begin // Leading 1 at bit 30
                leading_zeros = 5'd16;
                normalized_mantissa = mant_prod[29:7];
                final_exponent = p1_exp - 8'd16;
            end
            25'b0000000000000000001??????: begin // Leading 1 at bit 29
                leading_zeros = 5'd17;
                normalized_mantissa = mant_prod[28:6];
                final_exponent = p1_exp - 8'd17;
            end
            25'b00000000000000000001?????: begin // Leading 1 at bit 28
                leading_zeros = 5'd18;
                normalized_mantissa = mant_prod[27:5];
                final_exponent = p1_exp - 8'd18;
            end
            25'b000000000000000000001????: begin // Leading 1 at bit 27
                leading_zeros = 5'd19;
                normalized_mantissa = mant_prod[26:4];
                final_exponent = p1_exp - 8'd19;
            end
            25'b0000000000000000000001???: begin // Leading 1 at bit 26
                leading_zeros = 5'd20;
                normalized_mantissa = mant_prod[25:3];
                final_exponent = p1_exp - 8'd20;
            end
            25'b00000000000000000000001??: begin // Leading 1 at bit 25
                leading_zeros = 5'd21;
                normalized_mantissa = mant_prod[24:2];
                final_exponent = p1_exp - 8'd21;
            end
            25'b000000000000000000000001?: begin // Leading 1 at bit 24
                leading_zeros = 5'd22;
                normalized_mantissa = mant_prod[23:1];
                final_exponent = p1_exp - 8'd22;
            end
            25'b0000000000000000000000001: begin // Leading 1 at bit 23
                leading_zeros = 5'd23;
                normalized_mantissa = mant_prod[22:0];
                final_exponent = p1_exp - 8'd23;
            end
            default: begin // All zeros or underflow
                leading_zeros = 5'd24;
                normalized_mantissa = 23'b0;
                final_exponent = 8'b0; // Zero result
            end
        endcase
    end
    
    // --- PIPELINE STAGE 2: Normalize and Assemble Final Result ---
    // This final stage takes the product and combines it with the
    // correctly synchronized sign and exponent from the previous stage's registers.
    always @(posedge clk or posedge rst) begin
        if (rst) begin
            result <= 32'b0;
        end else begin
            // Check for zero result or underflow
            if (final_exponent == 8'b0 || final_exponent[7] == 1'b1) begin
                // Zero result or exponent underflow
                result <= {p1_sign, 31'b0};
            end else if (final_exponent == 8'hFF) begin
                // Exponent overflow - infinity
                result <= {p1_sign, 8'hFF, 23'b0};
            end else begin
                // Normal result
                result <= {p1_sign, final_exponent, normalized_mantissa};
            end
        end
    end

endmodule
// Wallace Tree Multiplier for 24x24 bits
module mult(
    input [23:0] a,
    input [23:0] b,
    output [47:0] product
);

    // Partial products - 24x24 = 576 partial products
    wire [47:0] pp [0:23]; // 24 partial products, each 48 bits wide
    
    // Generate partial products
    genvar i, j;
    generate
        for (i = 0; i < 24; i = i + 1) begin : pp_gen
            for (j = 0; j < 24; j = j + 1) begin : bit_gen
                if (j + i < 48) begin
                    assign pp[i][j + i] = a[j] & b[i];
                end
            end
            // Fill unused bits with zeros
            for (j = 0; j < i; j = j + 1) begin : zero_fill_low
                assign pp[i][j] = 1'b0;
            end
            for (j = i + 24; j < 48; j = j + 1) begin : zero_fill_high
                assign pp[i][j] = 1'b0;
            end
        end
    endgenerate

    // Wallace Tree Reduction
    // Level 1: Reduce 24 partial products to 16
    wire [47:0] level1 [0:15];
    
    // Full adders for 3:2 compression
    generate
        for (i = 0; i < 16; i = i + 1) begin : level1_reduction
            if (i < 8) begin
                // Use 3:2 compressors (full adders)
                wallace_csa_48bit csa1_l1 (
                    .a(pp[i*3]),
                    .b(pp[i*3 + 1]),
                    .c(pp[i*3 + 2]),
                    .sum(level1[i*2]),
                    .carry(level1[i*2 + 1])
                );
            end
        end
    endgenerate
    
    // Level 2: Reduce 16 to 11
    wire [47:0] level2 [0:10];
    
    generate
        for (i = 0; i < 5; i = i + 1) begin : level2_reduction
            wallace_csa_48bit csa2_l2 (
                .a(level1[i*3]),
                .b(level1[i*3 + 1]),
                .c(level1[i*3 + 2]),
                .sum(level2[i*2]),
                .carry(level2[i*2 + 1])
            );
        end
        // Handle remaining partial product
        assign level2[10] = level1[15];
    endgenerate
    
    // Level 3: Reduce 11 to 8
    wire [47:0] level3 [0:7];
    
    generate
        for (i = 0; i < 3; i = i + 1) begin : level3_reduction
            wallace_csa_48bit csa3_l3 (
                .a(level2[i*3]),
                .b(level2[i*3 + 1]),
                .c(level2[i*3 + 2]),
                .sum(level3[i*2]),
                .carry(level3[i*2 + 1])
            );
        end
        // Handle remaining partial products
        assign level3[6] = level2[9];
        assign level3[7] = level2[10];
    endgenerate
    
    // Level 4: Reduce 8 to 6
    wire [47:0] level4 [0:5];
    
    generate
        for (i = 0; i < 2; i = i + 1) begin : level4_reduction
            wallace_csa_48bit csa4_l4 (
                .a(level3[i*3]),
                .b(level3[i*3 + 1]),
                .c(level3[i*3 + 2]),
                .sum(level4[i*2]),
                .carry(level4[i*2 + 1])
            );
        end
        // Handle remaining partial products
        assign level4[4] = level3[6];
        assign level4[5] = level3[7];
    endgenerate
    
    // Level 5: Reduce 6 to 4
    wire [47:0] level5 [0:3];
    
    generate
        for (i = 0; i < 2; i = i + 1) begin : level5_reduction
            wallace_csa_48bit csa5_l5 (
                .a(level4[i*3]),
                .b(level4[i*3 + 1]),
                .c(level4[i*3 + 2]),
                .sum(level5[i*2]),
                .carry(level5[i*2 + 1])
            );
        end
    endgenerate
    
    // Level 6: Reduce 4 to 3
    wire [47:0] level6 [0:2];
    
    wallace_csa_48bit csa6_l6 (
        .a(level5[0]),
        .b(level5[1]),
        .c(level5[2]),
        .sum(level6[0]),
        .carry(level6[1])
    );
    assign level6[2] = level5[3];
    
    // Level 7: Reduce 3 to 2
    wire [47:0] final_sum, final_carry;
    
    wallace_csa_48bit csa7_final (
        .a(level6[0]),
        .b(level6[1]),
        .c(level6[2]),
        .sum(final_sum),
        .carry(final_carry)
    );
    
    // Final addition using Kogge-Stone adder
    wire [47:0] shifted_carry;
    assign shifted_carry = {final_carry[46:0], 1'b0}; // Left shift carry by 1
    
    // Use multiple 8-bit Kogge-Stone adders for final 48-bit addition
    wire [5:0] carry_chain; // Carry between 8-bit adders
    
    kogge_stone_48bit final_adder (
        .a(final_sum),
        .b(shifted_carry),
        .sum(product)
    );

endmodule

// 48-bit Carry-Save Adder (3:2 compressor)
module wallace_csa_48bit(
    input [47:0] a,
    input [47:0] b,
    input [47:0] c,
    output [47:0] sum,
    output [47:0] carry
);

    genvar i;
    generate
        for (i = 0; i < 48; i = i + 1) begin : csa_bits
            full_adder fa (
                .a(a[i]),
                .b(b[i]),
                .cin(c[i]),
                .sum(sum[i]),
                .cout(carry[i])
            );
        end
    endgenerate

endmodule

// 48-bit Kogge-Stone Adder
module kogge_stone_48bit(
    input [47:0] a,
    input [47:0] b,
    output [47:0] sum
);

    // Use six 8-bit Kogge-Stone adders with carry propagation
    wire [5:0] carry_out;
    
    kogge_stone8 add0 (.x(a[7:0]),   .y(b[7:0]),   .sum(sum[7:0]),   .cin(1'b0),        .cout(carry_out[0]));
    kogge_stone8 add1 (.x(a[15:8]),  .y(b[15:8]),  .sum(sum[15:8]),  .cin(carry_out[0]), .cout(carry_out[1]));
    kogge_stone8 add2 (.x(a[23:16]), .y(b[23:16]), .sum(sum[23:16]), .cin(carry_out[1]), .cout(carry_out[2]));
    kogge_stone8 add3 (.x(a[31:24]), .y(b[31:24]), .sum(sum[31:24]), .cin(carry_out[2]), .cout(carry_out[3]));
    kogge_stone8 add4 (.x(a[39:32]), .y(b[39:32]), .sum(sum[39:32]), .cin(carry_out[3]), .cout(carry_out[4]));
    kogge_stone8 add5 (.x(a[47:40]), .y(b[47:40]), .sum(sum[47:40]), .cin(carry_out[4]), .cout(carry_out[5]));

endmodule

// Full Adder
module full_adder(
    input a,
    input b,
    input cin,
    output sum,
    output cout
);
    assign sum = a ^ b ^ cin;
    assign cout = (a & b) | (cin & (a ^ b));
endmodule


module kogge_stone8(x, y, sum, cin, cout);
    input [7:0] y, x;
    output [7:0] sum;
    input cin;
    output cout;

    wire [7:0] G_Z, P_Z;
    wire [7:0] G_A, P_A;
    wire [7:0] G_B, P_B;
    wire [7:0] G_C, P_C;

    // level 0
    genvar i;
    generate
        for (i = 0; i < 8; i = i + 1) begin : gen0
            and_xor ax(x[i], y[i], P_Z[i], G_Z[i]);
        end
    endgenerate

    // level 1
    gray_cell gc0(cin, P_Z[0], G_Z[0], G_A[0]);
    black_cell bc1(G_Z[0], P_Z[1], G_Z[1], P_Z[0], G_A[1], P_A[1]);
    black_cell bc2(G_Z[1], P_Z[2], G_Z[2], P_Z[1], G_A[2], P_A[2]);
    black_cell bc3(G_Z[2], P_Z[3], G_Z[3], P_Z[2], G_A[3], P_A[3]);
    black_cell bc4(G_Z[3], P_Z[4], G_Z[4], P_Z[3], G_A[4], P_A[4]);
    black_cell bc5(G_Z[4], P_Z[5], G_Z[5], P_Z[4], G_A[5], P_A[5]);
    black_cell bc6(G_Z[5], P_Z[6], G_Z[6], P_Z[5], G_A[6], P_A[6]);
    black_cell bc7(G_Z[6], P_Z[7], G_Z[7], P_Z[6], G_A[7], P_A[7]);

    // level 2
    gray_cell gc1(cin, P_A[1], G_A[1], G_B[1]);
    gray_cell gc2(G_A[0], P_A[2], G_A[2], G_B[2]);
    black_cell bc3b(G_A[1], P_A[3], G_A[3], P_A[1], G_B[3], P_B[3]);
    black_cell bc4b(G_A[2], P_A[4], G_A[4], P_A[2], G_B[4], P_B[4]);
    black_cell bc5b(G_A[3], P_A[5], G_A[5], P_A[3], G_B[5], P_B[5]);
    black_cell bc6b(G_A[4], P_A[6], G_A[6], P_A[4], G_B[6], P_B[6]);
    black_cell bc7b(G_A[5], P_A[7], G_A[7], P_A[5], G_B[7], P_B[7]);

    // level 3
    gray_cell gc3(cin, P_B[3], G_B[3], G_C[3]);
    gray_cell gc4(G_A[0], P_B[4], G_B[4], G_C[4]);
    gray_cell gc5(G_B[1], P_B[5], G_B[5], G_C[5]);
    gray_cell gc6(G_B[2], P_B[6], G_B[6], G_C[6]);
    black_cell bc7c(G_B[3], P_B[7], G_B[7], P_B[3], G_C[7], P_C[7]);

    // level 4
    gray_cell gc7(cin, P_C[7], G_C[7], cout);

    // sum generation
    xor(sum[0], cin, P_Z[0]);
    xor(sum[1], G_A[0], P_Z[1]);
    xor(sum[2], G_B[1], P_Z[2]);
    xor(sum[3], G_B[2], P_Z[3]);
    xor(sum[4], G_C[3], P_Z[4]);
    xor(sum[5], G_C[4], P_Z[5]);
    xor(sum[6], G_C[5], P_Z[6]);
    xor(sum[7], G_C[6], P_Z[7]);
endmodule

//now we have cmp, addsub, muldiv take a control signal, clk rst ,a,b and consig must be 3 bit long, compare, add, sub, divide

module fpu_top(
    input        clk,
    input        rst,
    input [31:0] a,
    input [31:0] b,
    input  [2:0] consig, // 3-bit control: 000=cmp, 001=add, 010=sub, 011=mul, 100=div (not implemented)
    output reg [31:0] out,
    output reg [2:0] cmp_out
);

    wire [31:0] addsub_out;
    wire [2:0]  cmp_con;
    wire [31:0] muldiv_out;

    // Instantiate compare
    cmp cmp_inst (
        .clk(clk),
        .rst(rst),
        .a(a),
        .b(b),
        .con(cmp_con)
    );

    // Instantiate add/sub
    addsub addsub_inst (
        .clk(clk),
        .rst(rst),
        .a(a),
        .b(b),
        .x(consig[0]), // 0: add, 1: sub
        .out(addsub_out)
    );

    // Instantiate mul
    muldiv muldiv_inst (
        .clk(clk),
        .rst(rst),
        .a(a),
        .b(b),
        .result(muldiv_out)
    );

    always @(*) begin
        case (consig)
            3'b000: begin // compare
                out = 32'b0;
                cmp_out = cmp_con;
            end
            3'b001: begin // add
                out = addsub_out;
                cmp_out = 3'b0;
            end
            3'b010: begin // sub
                out = addsub_out;
                cmp_out = 3'b0;
            end
            3'b011: begin // mul
                out = muldiv_out;
                cmp_out = 3'b0;
            end
            default: begin
                out = 32'b0;
                cmp_out = 3'b0;
            end
        endcase
    end

endmodule

