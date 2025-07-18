module fpu_tb;

    // Testbench signals
    reg clk;
    reg rst;
    reg [31:0] a;
    reg [31:0] b;
    reg [2:0] consig;
    wire [31:0] out;
    wire [2:0] cmp_out;
    
    // Test vectors storage
    reg [31:0] test_a [0:15];
    reg [31:0] test_b [0:15];
    reg [2:0] test_ops [0:4];
    
    // Loop variables
    integer i, j;
    
    // Clock generation
    initial begin
        clk = 0;
        forever #5 clk = ~clk; // 100MHz clock (10ns period)
    end
    
    // Instantiate the FPU
    fpu_top dut (
        .clk(clk),
        .rst(rst),
        .a(a),
        .b(b),
        .consig(consig),
        .out(out),
        .cmp_out(cmp_out)
    );
    
    // Initialize test vectors
    initial begin
        // IEEE 754 test values
        test_a[0]  = 32'h3F800000;  // 1.0
        test_a[1]  = 32'h40000000;  // 2.0
        test_a[2]  = 32'h40400000;  // 3.0
        test_a[3]  = 32'h40800000;  // 4.0
        test_a[4]  = 32'h40A00000;  // 5.0
        test_a[5]  = 32'hBF800000;  // -1.0
        test_a[6]  = 32'hC0000000;  // -2.0
        test_a[7]  = 32'h41200000;  // 10.0
        test_a[8]  = 32'h3F000000;  // 0.5
        test_a[9]  = 32'h3F400000;  // 0.75
        test_a[10] = 32'h00000000;  // 0.0
        test_a[11] = 32'h80000000;  // -0.0
        test_a[12] = 32'h42C80000;  // 100.0
        test_a[13] = 32'h3DCCCCCD;  // 0.1
        test_a[14] = 32'h7F7FFFFF;  // Max normal
        test_a[15] = 32'h00800000;  // Min normal
        
        test_b[0]  = 32'h40000000;  // 2.0
        test_b[1]  = 32'h3F800000;  // 1.0
        test_b[2]  = 32'h40000000;  // 2.0
        test_b[3]  = 32'h40000000;  // 2.0
        test_b[4]  = 32'h40400000;  // 3.0
        test_b[5]  = 32'h40000000;  // 2.0
        test_b[6]  = 32'h3F800000;  // 1.0
        test_b[7]  = 32'h40000000;  // 2.0
        test_b[8]  = 32'h40000000;  // 2.0
        test_b[9]  = 32'h3F800000;  // 1.0
        test_b[10] = 32'h3F800000;  // 1.0
        test_b[11] = 32'h00000000;  // 0.0
        test_b[12] = 32'h41200000;  // 10.0
        test_b[13] = 32'h3F800000;  // 1.0
        test_b[14] = 32'h3F800000;  // 1.0
        test_b[15] = 32'h40000000;  // 2.0
        
        // Operations: compare, add, sub, mul, div
        test_ops[0] = 3'b000; // Compare
        test_ops[1] = 3'b001; // Add
        test_ops[2] = 3'b010; // Subtract
        test_ops[3] = 3'b011; // Multiply
        test_ops[4] = 3'b100; // Divide (if implemented)
    end
    
    // Helper function to convert IEEE 754 to real (for display)
    function real ieee754_to_real;
        input [31:0] ieee;
        reg sign;
        reg [7:0] exponent;
        reg [22:0] mantissa;
        real result;
        begin
            sign = ieee[31];
            exponent = ieee[30:23];
            mantissa = ieee[22:0];
            
            if (exponent == 8'h00) begin
                if (mantissa == 23'h0) begin
                    result = 0.0; // Zero
                end else begin
                    result = 0.0; // Denormalized (simplified)
                end
            end else if (exponent == 8'hFF) begin
                result = 999999.0; // Infinity/NaN (simplified)
            end else begin
                // Normal number
                result = (1.0 + mantissa / (2.0 ** 23)) * (2.0 ** (exponent - 127));
                if (sign) result = -result;
            end
            ieee754_to_real = result;
        end
    endfunction
    
    // Helper function to get operation name
    function [63:0] get_op_name;
        input [2:0] op;
        begin
            case (op)
                3'b000: get_op_name = "COMPARE ";
                3'b001: get_op_name = "ADD     ";
                3'b010: get_op_name = "SUBTRACT";
                3'b011: get_op_name = "MULTIPLY";
                3'b100: get_op_name = "DIVIDE  ";
                default: get_op_name = "UNKNOWN ";
            endcase
        end
    endfunction
    
    // Main test sequence
    initial begin
        // Initialize
        $display("=====================================");
        $display("       FPU Testbench Started");
        $display("=====================================");
        
        rst = 1;
        a = 0;
        b = 0;
        consig = 0;
        
        // Wait for reset
        #20;
        rst = 0;
        #10;
        
        // Test each operation with different inputs
        for (j = 0; j < 5; j = j + 1) begin
            $display("\n--- Testing %s ---", get_op_name(test_ops[j]));
            consig = test_ops[j];
            
            for (i = 0; i < 8; i = i + 1) begin
                a = test_a[i];
                b = test_b[i];
                
                // Wait for computation (pipelined operations need multiple cycles)
                #30;
                
                $display("Test %0d: A=%h (%.3f), B=%h (%.3f)", 
                    i+1, a, ieee754_to_real(a), b, ieee754_to_real(b));
                
                if (consig == 3'b000) begin
                    // Compare operation
                    $display("  Result: cmp_out=%b", cmp_out);
                    $display("  Signs: A_sign=%b, B_sign=%b, Equal_mag=%b", 
                        cmp_out[2], cmp_out[1], cmp_out[0]);
                end else begin
                    // Arithmetic operations
                    $display("  Result: %h (%.3f)", out, ieee754_to_real(out));
                end
                
                $display("  --------------------------------");
            end
        end
        
        // Specific test cases
        $display("\n=== Specific Test Cases ===");
        
        // Test 1: Simple addition (1.0 + 2.0 = 3.0)
        $display("\nTest: 1.0 + 2.0 = ?");
        consig = 3'b001; // Add
        a = 32'h3F800000; // 1.0
        b = 32'h40000000; // 2.0
        #30;
        $display("Expected: 3.0 (40400000), Got: %.3f (%h)", ieee754_to_real(out), out);
        
        // Test 2: Simple subtraction (3.0 - 1.0 = 2.0)
        $display("\nTest: 3.0 - 1.0 = ?");
        consig = 3'b010; // Subtract
        a = 32'h40400000; // 3.0
        b = 32'h3F800000; // 1.0
        #30;
        $display("Expected: 2.0 (40000000), Got: %.3f (%h)", ieee754_to_real(out), out);
        
        // Test 3: Simple multiplication (2.0 * 3.0 = 6.0)
        $display("\nTest: 2.0 * 3.0 = ?");
        consig = 3'b011; // Multiply
        a = 32'h40000000; // 2.0
        b = 32'h40400000; // 3.0
        #30;
        $display("Expected: 6.0 (40C00000), Got: %.3f (%h)", ieee754_to_real(out), out);
        
        // Test 4: Zero handling
        $display("\nTest: 0.0 + 5.0 = ?");
        consig = 3'b001; // Add
        a = 32'h00000000; // 0.0
        b = 32'h40A00000; // 5.0
        #30;
        $display("Expected: 5.0 (40A00000), Got: %.3f (%h)", ieee754_to_real(out), out);
        
        // Test 5: Negative numbers
        $display("\nTest: -1.0 + 2.0 = ?");
        consig = 3'b001; // Add
        a = 32'hBF800000; // -1.0
        b = 32'h40000000; // 2.0
        #30;
        $display("Expected: 1.0 (3F800000), Got: %.3f (%h)", ieee754_to_real(out), out);
        
        // Test 6: Compare equal numbers
        $display("\nTest: Compare 2.0 vs 2.0");
        consig = 3'b000; // Compare
        a = 32'h40000000; // 2.0
        b = 32'h40000000; // 2.0
        #30;
        $display("Expected: Equal (cmp_out[0]=1), Got: cmp_out=%b", cmp_out);
        
        // Test 7: Compare different signs
        $display("\nTest: Compare 2.0 vs -2.0");
        consig = 3'b000; // Compare
        a = 32'h40000000; // 2.0
        b = 32'hC0000000; // -2.0
        #30;
        $display("Got: cmp_out=%b (A_sign=%b, B_sign=%b)", 
            cmp_out, cmp_out[2], cmp_out[1]);
        
        // End simulation
        #50;
        $display("\n=====================================");
        $display("       FPU Testbench Completed");
        $display("=====================================");
        $finish;
    end
    
    // Monitor signals for debugging
    initial begin
        $monitor("Time=%0t, consig=%b, a=%h, b=%h, out=%h, cmp_out=%b", 
                 $time, consig, a, b, out, cmp_out);
    end
    
    // Generate VCD file for waveform viewing
    initial begin
        $dumpfile("fpu_tb.vcd");
        $dumpvars(0, fpu_tb);
    end
    
endmodule

