// terms: 7
// fval:  79374.25

module unsigned_exchange_8x8_l4_lamb8000_7 (
	input [7:0] x,
	input [7:0] y,
	output [15:0] z
);

wire [7:0] part1 =  y & {8{x[0]}};
wire [7:0] part2 =  y & {8{x[1]}};
wire [7:0] part3 =  y & {8{x[2]}};
wire [7:0] part4 =  y & {8{x[3]}};
wire [7:0] part5 =  y & {8{x[4]}};
wire [7:0] part6 =  y & {8{x[5]}};
wire [7:0] part7 =  y & {8{x[6]}};
wire [7:0] part8 =  y & {8{x[7]}};

wire [10:0] new_part1;
assign new_part1[0] = 0;
assign new_part1[1] = 0;
assign new_part1[2] = 0;
assign new_part1[3] = 0;
assign new_part1[4] = 0;
assign new_part1[5] = 0;
assign new_part1[6] = 0;
assign new_part1[7] = part1[7] | part2[6];
assign new_part1[8] = part2[7];
assign new_part1[9] = part3[6] | part4[5];
assign new_part1[10] = part3[7] & part4[6];

wire [10:0] new_part2;
assign new_part2[0] = 0;
assign new_part2[1] = 0;
assign new_part2[2] = 0;
assign new_part2[3] = 0;
assign new_part2[4] = 0;
assign new_part2[5] = 0;
assign new_part2[6] = 0;
assign new_part2[7] = part3[5] | part4[4];
assign new_part2[8] = 0;
assign new_part2[9] = part3[7] ^ part4[6];
assign new_part2[10] = part4[7];

wire [11:0] tmp_z = y*x[7:4];

assign z = {tmp_z, 4'd 0} + new_part1 + new_part2;

endmodule
