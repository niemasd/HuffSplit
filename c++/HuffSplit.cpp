/* AUTHOR: Niema Moshiri
 * DNA Split Huffman Compression
 *
 * USAGE:
 * -Compress:   HuffSplit compress <in_file>
 * -Decompress: HuffSplit decompress <huffsplit_file>
 *
 * COMPRESSED FILE OUTPUT FORMAT:
 * -The first byte of the compressed file ("InfoByte") tells us the tree topology (out of the 165 possible topologies)
 * -The next 4 bytes of the compressed file represent an int telling us how many total symbols are in the compressed file
 * -The remaining bytes are the compressed data
 *
 * If there is only 1 unique symbol, the resulting compressed file will only contain the first 9 bytes ("InfoByte" + "numChars")
 *
 * NOTE: Because of the DP algorithm to find optimal cuts, the message MUST be able to fit comfortably into RAM!!!
 */
#include<iostream>
using namespace std;
int main() {
    cout << "HI" << endl;
}