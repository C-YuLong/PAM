// terms: 3
// fval:  8939.25

#include <bitset>
#include <cmath>

using namespace std;

uint16_t unsigned_exchange_8x8_l2_lamb2000_0(uint16_t i, uint16_t j) {

bitset<8> x = i;

bitset<8> part1 = j*x[0];
bitset<8> part2 = j*x[1];
bitset<8> part3 = j*x[2];
bitset<8> part4 = j*x[3];
bitset<8> part5 = j*x[4];
bitset<8> part6 = j*x[5];
bitset<8> part7 = j*x[6];
bitset<8> part8 = j*x[7];

uint16_t tmp_z =  j*x[2]*pow(2, 2) + j*x[3]*pow(2, 3) + j*x[4]*pow(2, 4) + j*x[5]*pow(2, 5) + j*x[6]*pow(2, 6) + j*x[7]*pow(2, 7);

uint16_t z = tmp_z + (part1[6]|part2[5])*pow(2, 7) + (part1[7]|part2[6])*pow(2, 7) + (part2[7]*pow(2, 8));

return z;

}
