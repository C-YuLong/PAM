// terms: 7
// fval:  85394.25

#include <bitset>
#include <cmath>

using namespace std;

uint16_t unsigned_exchange_8x8_l4_lamb9000_9(uint16_t i, uint16_t j) {

bitset<8> x = i;

bitset<8> part1 = j*x[0];
bitset<8> part2 = j*x[1];
bitset<8> part3 = j*x[2];
bitset<8> part4 = j*x[3];
bitset<8> part5 = j*x[4];
bitset<8> part6 = j*x[5];
bitset<8> part7 = j*x[6];
bitset<8> part8 = j*x[7];

uint16_t tmp_z =  j*x[4]*pow(2, 4) + j*x[5]*pow(2, 5) + j*x[6]*pow(2, 6) + j*x[7]*pow(2, 7);

uint16_t z = tmp_z + (part1[7]|part2[6])*pow(2, 8) + (part2[7]*pow(2, 8)) + (part3[6]|part4[4])*pow(2, 8) + (part3[5]|part4[5])*pow(2, 8) + (part3[7]&part4[6])*pow(2, 10) + (part3[7]|part4[6])*pow(2, 9) + (part4[7]*pow(2, 10));

return z;

}
