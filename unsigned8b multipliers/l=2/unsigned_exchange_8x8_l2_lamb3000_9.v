// terms: 4
// fval:  14607.25

module unsigned_exchange_8x8_l2_lamb3000_9 (
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

wire [8:0] new_part1;
assign new_part1[0] = 0;
assign new_part1[1] = 0;
assign new_part1[2] = 0;
assign new_part1[3] = 0;
assign new_part1[4] = 0;
assign new_part1[5] = 0;
assign new_part1[6] = part1[5] | part2[4];
assign new_part1[7] = part1[7] | part2[6];
assign new_part1[8] = part1[7] & part2[6];

wire [8:0] new_part2;
assign new_part2[0] = 0;
assign new_part2[1] = 0;
assign new_part2[2] = 0;
assign new_part2[3] = 0;
assign new_part2[4] = 0;
assign new_part2[5] = 0;
assign new_part2[6] = 0;
assign new_part2[7] = 0;
assign new_part2[8] = part2[7];

wire [13:0] tmp_z = y*x[7:2];

assign z = {tmp_z, 2'd 0} + new_part1 + new_part2;

endmodule
