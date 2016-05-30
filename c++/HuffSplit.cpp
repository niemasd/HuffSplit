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
// definitions
#ifndef NUMTOPS
#define NUMTOPS 165
#endif
typedef unsigned char byte;

// includes
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <stdlib.h>
using namespace std;

// global variables
vector<unordered_map<char,string>> TOPS (NUMTOPS);

/* Helper Class: Node
 */
class Node {
    public:
        char symbol;
        Node* c0 = NULL;
        Node* c1 = NULL;
        Node( char s );
};
Node::Node( char s ) {
    symbol = s;
}

/* Given an integer, return the corresponding topology (see Topologies.pptx)
 * INPUT:  An integer (between 0 and 164, inclusive)
 * OUTPUT: A HashMap containing the code of the corresponding topology
 */
unordered_map<char,string> getCode( int topology ) {
    unordered_map<char,string> code;
    switch(topology) {
        // 1 Unique Character
        case 0:   code['A'] = ""; break;
        case 1:   code['C'] = ""; break;
        case 2:   code['G'] = ""; break;
        case 3:   code['T'] = ""; break;
        case 4:   code['N'] = ""; break;
        
        // 2 Unique Characters
        case 5:   code['A'] = "1"; code['C'] = "0"; break;
        case 6:   code['A'] = "1"; code['G'] = "0"; break;
        case 7:   code['A'] = "1"; code['T'] = "0"; break;
        case 8:   code['A'] = "1"; code['N'] = "0"; break;
        case 9:   code['C'] = "1"; code['G'] = "0"; break;
        case 10:  code['C'] = "1"; code['T'] = "0"; break;
        case 11:  code['C'] = "1"; code['N'] = "0"; break;
        case 12:  code['G'] = "1"; code['T'] = "0"; break;
        case 13:  code['G'] = "1"; code['N'] = "0"; break;
        case 14:  code['T'] = "1"; code['N'] = "0"; break;
        
        // 3 Unique Characters
        case 15:  code['A'] = "1"; code['C'] = "01"; code['G'] = "00"; break;
        case 16:  code['A'] = "1"; code['C'] = "01"; code['T'] = "00"; break;
        case 17:  code['A'] = "1"; code['C'] = "01"; code['N'] = "00"; break;
        case 18:  code['A'] = "1"; code['G'] = "01"; code['T'] = "00"; break;
        case 19:  code['A'] = "1"; code['G'] = "01"; code['N'] = "00"; break;
        case 20:  code['A'] = "1"; code['T'] = "01"; code['N'] = "00"; break;
        case 21:  code['C'] = "1"; code['A'] = "01"; code['G'] = "00"; break;
        case 22:  code['C'] = "1"; code['A'] = "01"; code['T'] = "00"; break;
        case 23:  code['C'] = "1"; code['A'] = "01"; code['N'] = "00"; break;
        case 24:  code['C'] = "1"; code['G'] = "01"; code['T'] = "00"; break;
        case 25:  code['C'] = "1"; code['G'] = "01"; code['N'] = "00"; break;
        case 26:  code['C'] = "1"; code['T'] = "01"; code['N'] = "00"; break;
        case 27:  code['G'] = "1"; code['A'] = "01"; code['C'] = "00"; break;
        case 28:  code['G'] = "1"; code['A'] = "01"; code['T'] = "00"; break;
        case 29:  code['G'] = "1"; code['A'] = "01"; code['N'] = "00"; break;
        case 30:  code['G'] = "1"; code['C'] = "01"; code['T'] = "00"; break;
        case 31:  code['G'] = "1"; code['C'] = "01"; code['N'] = "00"; break;
        case 32:  code['G'] = "1"; code['T'] = "01"; code['N'] = "00"; break;
        case 33:  code['T'] = "1"; code['A'] = "01"; code['C'] = "00"; break;
        case 34:  code['T'] = "1"; code['A'] = "01"; code['G'] = "00"; break;
        case 35:  code['T'] = "1"; code['A'] = "01"; code['N'] = "00"; break;
        case 36:  code['T'] = "1"; code['C'] = "01"; code['G'] = "00"; break;
        case 37:  code['T'] = "1"; code['C'] = "01"; code['N'] = "00"; break;
        case 38:  code['T'] = "1"; code['G'] = "01"; code['N'] = "00"; break;
        case 39:  code['N'] = "1"; code['A'] = "01"; code['C'] = "00"; break;
        case 40:  code['N'] = "1"; code['A'] = "01"; code['G'] = "00"; break;
        case 41:  code['N'] = "1"; code['A'] = "01"; code['T'] = "00"; break;
        case 42:  code['N'] = "1"; code['C'] = "01"; code['G'] = "00"; break;
        case 43:  code['N'] = "1"; code['C'] = "01"; code['T'] = "00"; break;
        case 44:  code['N'] = "1"; code['G'] = "01"; code['T'] = "00"; break;
        
        // 4 Unique Characters (Balanced)
        case 45:  code['A'] = "11"; code['C'] = "10"; code['G'] = "01"; code['T'] = "00"; break;
        case 46:  code['A'] = "11"; code['C'] = "10"; code['G'] = "01"; code['N'] = "00"; break;
        case 47:  code['A'] = "11"; code['C'] = "10"; code['T'] = "01"; code['N'] = "00"; break;
        case 48:  code['A'] = "11"; code['G'] = "10"; code['T'] = "01"; code['N'] = "00"; break;
        case 49:  code['C'] = "11"; code['G'] = "10"; code['T'] = "01"; code['N'] = "00"; break;
        
        // 4 Unique Characters (Unbalanced)
        case 50:  code['A'] = "1"; code['C'] = "01"; code['G'] = "001"; code['T'] = "000"; break;
        case 51:  code['A'] = "1"; code['C'] = "01"; code['G'] = "001"; code['N'] = "000"; break;
        case 52:  code['A'] = "1"; code['G'] = "01"; code['C'] = "001"; code['T'] = "000"; break;
        case 53:  code['A'] = "1"; code['G'] = "01"; code['C'] = "001"; code['N'] = "000"; break;
        case 54:  code['A'] = "1"; code['T'] = "01"; code['C'] = "001"; code['G'] = "000"; break;
        case 55:  code['A'] = "1"; code['T'] = "01"; code['C'] = "001"; code['N'] = "000"; break;
        case 56:  code['A'] = "1"; code['N'] = "01"; code['C'] = "001"; code['G'] = "000"; break;
        case 57:  code['A'] = "1"; code['N'] = "01"; code['C'] = "001"; code['T'] = "000"; break;
        case 58:  code['C'] = "1"; code['A'] = "01"; code['G'] = "001"; code['T'] = "000"; break;
        case 59:  code['C'] = "1"; code['A'] = "01"; code['G'] = "001"; code['N'] = "000"; break;
        case 60:  code['C'] = "1"; code['G'] = "01"; code['A'] = "001"; code['T'] = "000"; break;
        case 61:  code['C'] = "1"; code['G'] = "01"; code['A'] = "001"; code['N'] = "000"; break;
        case 62:  code['C'] = "1"; code['T'] = "01"; code['A'] = "001"; code['G'] = "000"; break;
        case 63:  code['C'] = "1"; code['T'] = "01"; code['A'] = "001"; code['N'] = "000"; break;
        case 64:  code['C'] = "1"; code['N'] = "01"; code['A'] = "001"; code['G'] = "000"; break;
        case 65:  code['C'] = "1"; code['N'] = "01"; code['A'] = "001"; code['T'] = "000"; break;
        case 66:  code['G'] = "1"; code['A'] = "01"; code['C'] = "001"; code['T'] = "000"; break;
        case 67:  code['G'] = "1"; code['A'] = "01"; code['C'] = "001"; code['N'] = "000"; break;
        case 68:  code['G'] = "1"; code['C'] = "01"; code['A'] = "001"; code['T'] = "000"; break;
        case 69:  code['G'] = "1"; code['C'] = "01"; code['A'] = "001"; code['N'] = "000"; break;
        case 70:  code['G'] = "1"; code['T'] = "01"; code['A'] = "001"; code['C'] = "000"; break;
        case 71:  code['G'] = "1"; code['T'] = "01"; code['A'] = "001"; code['N'] = "000"; break;
        case 72:  code['G'] = "1"; code['N'] = "01"; code['A'] = "001"; code['C'] = "000"; break;
        case 73:  code['G'] = "1"; code['N'] = "01"; code['A'] = "001"; code['T'] = "000"; break;
        case 74:  code['T'] = "1"; code['A'] = "01"; code['C'] = "001"; code['G'] = "000"; break;
        case 75:  code['T'] = "1"; code['A'] = "01"; code['C'] = "001"; code['N'] = "000"; break;
        case 76:  code['T'] = "1"; code['C'] = "01"; code['A'] = "001"; code['G'] = "000"; break;
        case 77:  code['T'] = "1"; code['C'] = "01"; code['A'] = "001"; code['N'] = "000"; break;
        case 78:  code['T'] = "1"; code['G'] = "01"; code['A'] = "001"; code['C'] = "000"; break;
        case 79:  code['T'] = "1"; code['G'] = "01"; code['A'] = "001"; code['N'] = "000"; break;
        case 80:  code['T'] = "1"; code['N'] = "01"; code['A'] = "001"; code['C'] = "000"; break;
        case 81:  code['T'] = "1"; code['N'] = "01"; code['A'] = "001"; code['G'] = "000"; break;
        case 82:  code['N'] = "1"; code['A'] = "01"; code['C'] = "001"; code['G'] = "000"; break;
        case 83:  code['N'] = "1"; code['A'] = "01"; code['C'] = "001"; code['T'] = "000"; break;
        case 84:  code['N'] = "1"; code['C'] = "01"; code['A'] = "001"; code['G'] = "000"; break;
        case 85:  code['N'] = "1"; code['C'] = "01"; code['A'] = "001"; code['T'] = "000"; break;
        case 86:  code['N'] = "1"; code['G'] = "01"; code['A'] = "001"; code['C'] = "000"; break;
        case 87:  code['N'] = "1"; code['G'] = "01"; code['A'] = "001"; code['T'] = "000"; break;
        case 88:  code['N'] = "1"; code['T'] = "01"; code['A'] = "001"; code['C'] = "000"; break;
        case 89:  code['N'] = "1"; code['T'] = "01"; code['A'] = "001"; code['G'] = "000"; break;
        
        // 5 Unique Characters (Line)
        case 90:  code['A'] = "1"; code['C'] = "01"; code['G'] = "001"; code['T'] = "0001"; code['N'] = "0000"; break;
        case 91:  code['A'] = "1"; code['C'] = "01"; code['T'] = "001"; code['G'] = "0001"; code['N'] = "0000"; break;
        case 92:  code['A'] = "1"; code['C'] = "01"; code['N'] = "001"; code['G'] = "0001"; code['T'] = "0000"; break;
        case 93:  code['A'] = "1"; code['G'] = "01"; code['C'] = "001"; code['T'] = "0001"; code['N'] = "0000"; break;
        case 94:  code['A'] = "1"; code['G'] = "01"; code['T'] = "001"; code['C'] = "0001"; code['N'] = "0000"; break;
        case 95:  code['A'] = "1"; code['G'] = "01"; code['N'] = "001"; code['C'] = "0001"; code['T'] = "0000"; break;
        case 96:  code['A'] = "1"; code['T'] = "01"; code['C'] = "001"; code['G'] = "0001"; code['N'] = "0000"; break;
        case 97:  code['A'] = "1"; code['T'] = "01"; code['G'] = "001"; code['C'] = "0001"; code['N'] = "0000"; break;
        case 98:  code['A'] = "1"; code['T'] = "01"; code['N'] = "001"; code['C'] = "0001"; code['G'] = "0000"; break;
        case 99:  code['A'] = "1"; code['N'] = "01"; code['C'] = "001"; code['G'] = "0001"; code['T'] = "0000"; break;
        case 100: code['A'] = "1"; code['N'] = "01"; code['G'] = "001"; code['C'] = "0001"; code['T'] = "0000"; break;
        case 101: code['A'] = "1"; code['N'] = "01"; code['T'] = "001"; code['C'] = "0001"; code['G'] = "0000"; break;
        case 102: code['C'] = "1"; code['A'] = "01"; code['G'] = "001"; code['T'] = "0001"; code['N'] = "0000"; break;
        case 103: code['C'] = "1"; code['A'] = "01"; code['T'] = "001"; code['G'] = "0001"; code['N'] = "0000"; break;
        case 104: code['C'] = "1"; code['A'] = "01"; code['N'] = "001"; code['G'] = "0001"; code['T'] = "0000"; break;
        case 105: code['C'] = "1"; code['G'] = "01"; code['A'] = "001"; code['T'] = "0001"; code['N'] = "0000"; break;
        case 106: code['C'] = "1"; code['G'] = "01"; code['T'] = "001"; code['A'] = "0001"; code['N'] = "0000"; break;
        case 107: code['C'] = "1"; code['G'] = "01"; code['N'] = "001"; code['A'] = "0001"; code['T'] = "0000"; break;
        case 108: code['C'] = "1"; code['T'] = "01"; code['A'] = "001"; code['G'] = "0001"; code['N'] = "0000"; break;
        case 109: code['C'] = "1"; code['T'] = "01"; code['G'] = "001"; code['A'] = "0001"; code['N'] = "0000"; break;
        case 110: code['C'] = "1"; code['T'] = "01"; code['N'] = "001"; code['A'] = "0001"; code['G'] = "0000"; break;
        case 111: code['C'] = "1"; code['N'] = "01"; code['A'] = "001"; code['G'] = "0001"; code['T'] = "0000"; break;
        case 112: code['C'] = "1"; code['N'] = "01"; code['G'] = "001"; code['A'] = "0001"; code['T'] = "0000"; break;
        case 113: code['C'] = "1"; code['N'] = "01"; code['T'] = "001"; code['A'] = "0001"; code['G'] = "0000"; break;
        case 114: code['G'] = "1"; code['A'] = "01"; code['C'] = "001"; code['T'] = "0001"; code['N'] = "0000"; break;
        case 115: code['G'] = "1"; code['A'] = "01"; code['T'] = "001"; code['C'] = "0001"; code['N'] = "0000"; break;
        case 116: code['G'] = "1"; code['A'] = "01"; code['N'] = "001"; code['C'] = "0001"; code['T'] = "0000"; break;
        case 117: code['G'] = "1"; code['C'] = "01"; code['A'] = "001"; code['T'] = "0001"; code['N'] = "0000"; break;
        case 118: code['G'] = "1"; code['C'] = "01"; code['T'] = "001"; code['A'] = "0001"; code['N'] = "0000"; break;
        case 119: code['G'] = "1"; code['C'] = "01"; code['N'] = "001"; code['A'] = "0001"; code['T'] = "0000"; break;
        case 120: code['G'] = "1"; code['T'] = "01"; code['A'] = "001"; code['C'] = "0001"; code['N'] = "0000"; break;
        case 121: code['G'] = "1"; code['T'] = "01"; code['C'] = "001"; code['A'] = "0001"; code['N'] = "0000"; break;
        case 122: code['G'] = "1"; code['T'] = "01"; code['N'] = "001"; code['A'] = "0001"; code['C'] = "0000"; break;
        case 123: code['G'] = "1"; code['N'] = "01"; code['A'] = "001"; code['C'] = "0001"; code['T'] = "0000"; break;
        case 124: code['G'] = "1"; code['N'] = "01"; code['C'] = "001"; code['A'] = "0001"; code['T'] = "0000"; break;
        case 125: code['G'] = "1"; code['N'] = "01"; code['T'] = "001"; code['A'] = "0001"; code['C'] = "0000"; break;
        case 126: code['T'] = "1"; code['A'] = "01"; code['C'] = "001"; code['G'] = "0001"; code['N'] = "0000"; break;
        case 127: code['T'] = "1"; code['A'] = "01"; code['G'] = "001"; code['C'] = "0001"; code['N'] = "0000"; break;
        case 128: code['T'] = "1"; code['A'] = "01"; code['N'] = "001"; code['C'] = "0001"; code['G'] = "0000"; break;
        case 129: code['T'] = "1"; code['C'] = "01"; code['A'] = "001"; code['G'] = "0001"; code['N'] = "0000"; break;
        case 130: code['T'] = "1"; code['C'] = "01"; code['G'] = "001"; code['A'] = "0001"; code['N'] = "0000"; break;
        case 131: code['T'] = "1"; code['C'] = "01"; code['N'] = "001"; code['A'] = "0001"; code['G'] = "0000"; break;
        case 132: code['T'] = "1"; code['G'] = "01"; code['A'] = "001"; code['C'] = "0001"; code['N'] = "0000"; break;
        case 133: code['T'] = "1"; code['G'] = "01"; code['C'] = "001"; code['A'] = "0001"; code['N'] = "0000"; break;
        case 134: code['T'] = "1"; code['G'] = "01"; code['N'] = "001"; code['A'] = "0001"; code['C'] = "0000"; break;
        case 135: code['T'] = "1"; code['N'] = "01"; code['A'] = "001"; code['C'] = "0001"; code['G'] = "0000"; break;
        case 136: code['T'] = "1"; code['N'] = "01"; code['C'] = "001"; code['A'] = "0001"; code['G'] = "0000"; break;
        case 137: code['T'] = "1"; code['N'] = "01"; code['G'] = "001"; code['A'] = "0001"; code['C'] = "0000"; break;
        case 138: code['N'] = "1"; code['A'] = "01"; code['C'] = "001"; code['G'] = "0001"; code['T'] = "0000"; break;
        case 139: code['N'] = "1"; code['A'] = "01"; code['G'] = "001"; code['C'] = "0001"; code['T'] = "0000"; break;
        case 140: code['N'] = "1"; code['A'] = "01"; code['T'] = "001"; code['C'] = "0001"; code['G'] = "0000"; break;
        case 141: code['N'] = "1"; code['C'] = "01"; code['A'] = "001"; code['G'] = "0001"; code['T'] = "0000"; break;
        case 142: code['N'] = "1"; code['C'] = "01"; code['G'] = "001"; code['A'] = "0001"; code['T'] = "0000"; break;
        case 143: code['N'] = "1"; code['C'] = "01"; code['T'] = "001"; code['A'] = "0001"; code['G'] = "0000"; break;
        case 144: code['N'] = "1"; code['G'] = "01"; code['A'] = "001"; code['C'] = "0001"; code['T'] = "0000"; break;
        case 145: code['N'] = "1"; code['G'] = "01"; code['C'] = "001"; code['A'] = "0001"; code['T'] = "0000"; break;
        case 146: code['N'] = "1"; code['G'] = "01"; code['T'] = "001"; code['A'] = "0001"; code['C'] = "0000"; break;
        case 147: code['N'] = "1"; code['T'] = "01"; code['A'] = "001"; code['C'] = "0001"; code['G'] = "0000"; break;
        case 148: code['N'] = "1"; code['T'] = "01"; code['C'] = "001"; code['A'] = "0001"; code['G'] = "0000"; break;
        case 149: code['N'] = "1"; code['T'] = "01"; code['G'] = "001"; code['A'] = "0001"; code['C'] = "0000"; break;
        
        // 5 Unique Characters (Bend 1)
        case 150: code['A'] = "11"; code['C'] = "10"; code['G'] = "01"; code['T'] = "001"; code['N'] = "000"; break;
        case 151: code['A'] = "11"; code['C'] = "10"; code['T'] = "01"; code['G'] = "001"; code['N'] = "000"; break;
        case 152: code['A'] = "11"; code['G'] = "10"; code['T'] = "01"; code['C'] = "001"; code['N'] = "000"; break;
        case 153: code['C'] = "11"; code['G'] = "10"; code['T'] = "01"; code['A'] = "001"; code['N'] = "000"; break;
        case 154: code['A'] = "11"; code['C'] = "10"; code['N'] = "01"; code['G'] = "001"; code['T'] = "000"; break;
        case 155: code['A'] = "11"; code['G'] = "10"; code['N'] = "01"; code['C'] = "001"; code['T'] = "000"; break;
        case 156: code['C'] = "11"; code['G'] = "10"; code['N'] = "01"; code['A'] = "001"; code['T'] = "000"; break;
        case 157: code['A'] = "11"; code['T'] = "10"; code['N'] = "01"; code['C'] = "001"; code['G'] = "000"; break;
        case 158: code['C'] = "11"; code['T'] = "10"; code['N'] = "01"; code['A'] = "001"; code['G'] = "000"; break;
        case 159: code['G'] = "11"; code['T'] = "10"; code['N'] = "01"; code['A'] = "001"; code['C'] = "000"; break;
        
        // 5 Unique Characters (Bend 2)
        case 160: code['A'] = "1"; code['C'] = "011"; code['G'] = "010"; code['T'] = "001"; code['N'] = "000"; break;
        case 161: code['C'] = "1"; code['A'] = "011"; code['G'] = "010"; code['T'] = "001"; code['N'] = "000"; break;
        case 162: code['G'] = "1"; code['A'] = "011"; code['C'] = "010"; code['T'] = "001"; code['N'] = "000"; break;
        case 163: code['T'] = "1"; code['A'] = "011"; code['C'] = "010"; code['G'] = "001"; code['N'] = "000"; break;
        case 164: code['N'] = "1"; code['A'] = "011"; code['C'] = "010"; code['G'] = "001"; code['T'] = "000"; break;
        
        // Out of Bounds
        default: cerr << "INVALID TREE TOPOLOGY: " << topology  << endl << endl; exit(-1);
    }
    return code;
}

/* Given an integer, return the corresponding tree (see Topologies.pptx)
 * INPUT:  An integer (between 0 and 164, inclusive)
 * OUTPUT: The root node of the corresponding tree
 */
Node* buildTree( int topology ) {
    unordered_map<char,string> code = getCode(topology);
    Node* root = new Node((char)0);
    if(topology < 5) {
        return root;
    }
    else if(topology < 15) {
        for(auto kv : code) {
            if(kv.second == "1") {
                root->c1 = new Node(kv.first);
            }
            else {
                root->c0 = new Node(kv.first);
            }
        }
    }
    else if(topology < 45) {
        Node* nodeA = NULL;
        Node* nodeB = NULL;
        Node* nodeC = NULL;
        for(auto kv : code) {
            if(kv.second == "1") {
                nodeA = new Node(kv.first);
            }
            else if(kv.second == "01") {
                nodeB = new Node(kv.first);
            }
            else {
                nodeC = new Node(kv.first);
            }
        }
        root->c1 = nodeA;
        root->c0 = new Node((char)0);
        root->c0->c1 = nodeB;
        root->c0->c0 = nodeC;
    }
    else if(topology < 50) {
        Node* nodeA = NULL;
        Node* nodeB = NULL;
        Node* nodeC = NULL;
        Node* nodeD = NULL;
        for(auto kv : code) {
            if(kv.second == "11") {
                nodeA = new Node(kv.first);
            }
            else if(kv.second == "10") {
                nodeB = new Node(kv.first);
            }
            else if(kv.second == "01") {
                nodeC = new Node(kv.first);
            }
            else {
                nodeD = new Node(kv.first);
            }
        }
        root->c1 = new Node((char)0);
        root->c0 = new Node((char)0);
        root->c1->c1 = nodeA;
        root->c1->c0 = nodeB;
        root->c0->c1 = nodeC;
        root->c0->c0 = nodeD;
    }
    else if(topology < 90) {
        Node* nodeA = NULL;
        Node* nodeB = NULL;
        Node* nodeC = NULL;
        Node* nodeD = NULL;
        for(auto kv : code) {
            if(kv.second == "1") {
                nodeA = new Node(kv.first);
            }
            else if(kv.second == "01") {
                nodeB = new Node(kv.first);
            }
            else if(kv.second == "001") {
                nodeC = new Node(kv.first);
            }
            else {
                nodeD = new Node(kv.first);
            }
        }
        root->c1 = nodeA;
        root->c0 = new Node((char)0);
        root->c0->c1 = nodeB;
        root->c0->c0 = new Node((char)0);
        root->c0->c0->c1 = nodeC;
        root->c0->c0->c0 = nodeD;
    }
    else if(topology < 150) {
        Node* nodeA = NULL;
        Node* nodeB = NULL;
        Node* nodeC = NULL;
        Node* nodeD = NULL;
        Node* nodeE = NULL;
        for(auto kv : code) {
            if(kv.second == "1") {
                nodeA = new Node(kv.first);
            }
            else if(kv.second == "01") {
                nodeB = new Node(kv.first);
            }
            else if(kv.second == "001") {
                nodeC = new Node(kv.first);
            }
            else if(kv.second == "0001") {
                nodeD = new Node(kv.first);
            }
            else {
                nodeE = new Node(kv.first);
            }
        }
        root->c1 = nodeA;
        root->c0 = new Node((char)0);
        root->c0->c1 = nodeB;
        root->c0->c0 = new Node((char)0);
        root->c0->c0->c1 = nodeC;
        root->c0->c0->c0 = new Node((char)0);
        root->c0->c0->c0->c1 = nodeD;
        root->c0->c0->c0->c0 = nodeE;
    }
    else if(topology < 160) {
        Node* nodeA = NULL;
        Node* nodeB = NULL;
        Node* nodeC = NULL;
        Node* nodeD = NULL;
        Node* nodeE = NULL;
        for(auto kv : code) {
            if(kv.second == "11") {
                nodeA = new Node(kv.first);
            }
            else if(kv.second == "10") {
                nodeB = new Node(kv.first);
            }
            else if(kv.second == "01") {
                nodeC = new Node(kv.first);
            }
            else if(kv.second == "001") {
                nodeD = new Node(kv.first);
            }
            else {
                nodeE = new Node(kv.first);
            }
        }
        root->c1 = new Node((char)0);
        root->c1->c1 = nodeA;
        root->c1->c0 = nodeB;
        root->c0 = new Node((char)0);
        root->c0->c1 = nodeC;
        root->c0->c0 = new Node((char)0);
        root->c0->c0->c1 = nodeD;
        root->c0->c0->c0 = nodeE;
    }
    else if(topology < 165) {
        Node* nodeA = NULL;
        Node* nodeB = NULL;
        Node* nodeC = NULL;
        Node* nodeD = NULL;
        Node* nodeE = NULL;
        for(auto kv : code) {
            if(kv.second == "1") {
                nodeA = new Node(kv.first);
            }
            else if(kv.second == "011") {
                nodeB = new Node(kv.first);
            }
            else if(kv.second == "010") {
                nodeC = new Node(kv.first);
            }
            else if(kv.second == "001") {
                nodeD = new Node(kv.first);
            }
            else {
                nodeE = new Node(kv.first);
            }
        }
        root->c1 = nodeA;
        root->c0 = new Node((char)0);
        root->c0->c1 = new Node((char)0);
        root->c0->c1->c1 = nodeB;
        root->c0->c1->c0 = nodeC;
        root->c0->c0 = new Node((char)0);
        root->c0->c0->c1 = nodeD;
        root->c0->c0->c0 = nodeE;
    }
    else {
        cerr << "ERROR: Trying to get tree for topology > 164!" << endl << endl; exit(-1);
    }
    return root;
}

/* Compress the input file using my split Huffman algorithm
 * INPUT:  A DNA string to compress
 * OUTPUT: The compressed results of my split Huffman algorithm
 */
int compress( string INFILE, string OUTFILE ) {
    // read input file as string 'in'
    ifstream infile(INFILE);
    string in(static_cast<stringstream const&>(stringstream() << infile.rdbuf()).str());
    char last = in[in.size()-1];
    if(last != 'A' && last != 'C' && last != 'G' && last != 'T' && last != 'N') {
        in = in.substr(0,in.size()-1);
    }
    const int L = in.size();

    // get optimal cuts
    int C[2][NUMTOPS];
    byte backtrack[L][NUMTOPS];
    int bestT[] = {-1,-1};
    int bestC[] = {-1,-1};
    for(int t = 0; t < NUMTOPS; ++t) {
        if(TOPS[t].find(in[0]) != TOPS[t].end()) {
            C[0][t] = 72 + TOPS[t][in[0]].size();
            backtrack[0][t] = t;
            if(bestC[0] == -1 || bestC[0] > C[0][t]) {
                bestC[0] = C[0][t];
                bestT[0] = t;
            }
        }
        else {
            C[0][t] = -1;
            backtrack[0][t] = 255;
        }
    }
    for(int i = 1; i < L; ++i) {
        int i1 = i%2;
        int i0 = (i-1)%2;
        char c = in[i];
        if(c != 'A' && c != 'C' && c != 'G' && c != 'T' && c != 'N') {
            cerr << "ERROR: Invalid symbol: " << c << endl << endl; exit(-1);
        }
        for(int top = 0; top < NUMTOPS; ++top) {
            if(TOPS[top].find(c) != TOPS[top].end()) {
                int bits = TOPS[top][c].size();
                if(bestT[i0] == top) {
                    C[i1][top] = C[i0][top] + bits;
                    backtrack[i][top] = top;
                }
                else {
                    int sameC = -1;
                    if(C[i0][top] != -1) {
                        sameC = C[i0][top] + bits;
                    }
                    int diffC = C[i0][bestT[i0]];
                    if(diffC%8 != 0) {
                        diffC += (8-(diffC%8));
                    }
                    diffC += (72+bits);
                    if(sameC == -1 || diffC < sameC) {
                        C[i1][top] = diffC;
                        backtrack[i][top] = bestT[i0];
                    }
                    else {
                        C[i1][top] = sameC;
                        backtrack[i][top] = top;
                    }
                }
                if(bestC[i1] == -1 || bestC[i1] > C[i1][top]) {
                    bestC[i1] = C[i1][top];
                    bestT[i1] = top;
                }
            }
            else {
                C[i1][top] = -1;
                backtrack[i][top] = 255;
            }
        }
        bestC[i0] = -1;
        bestT[i0] = -1;
    }

    // reconstruct topology path from backtrack
    byte path[L];
    const int LAST = (L-1)%2;
    for(int t = 0; t < NUMTOPS; ++t) {
        if(C[LAST][t] != -1 && (C[LAST][path[L-1]] == -1 || C[LAST][t] < C[LAST][path[L-1]])) {
            path[L-1] = t;
        }
    }
    for(int i = L-2; i >= 0; --i) {
        path[i] = backtrack[i+1][path[i+1]];
    }
    vector<int> cuts;
    cuts.push_back(0);
    for(int i = 1; i < L; ++i) {
        if(path[i] != path[i-1]) {
            cuts.push_back(i);
        }
    }
    cuts.push_back(L);

    // encode file
    ofstream out(OUTFILE,ios::binary);
    for(int cut = 1; cut < cuts.size(); ++cut) {
        // write header
        int start = cuts[cut-1];
        int end = cuts[cut];
        int length = end-start;
        out.write((char*)&path[start],sizeof(path[start]));
        out.write((char*)&length,sizeof(length));

        // if only 1 unique symbol, only need first 5 bytes
        if(path[start] < 5) {
            continue;
        }

        // encode substring
        string buf = "";
        for(int i = start; i < end; ++i) {
            buf += TOPS[path[start]][in[i]];
            if(buf.size() >= 8) {
                byte temp = (byte)strtoul(buf.c_str(),NULL,2);
                out.write((char*)&temp,sizeof(temp));
                buf = buf.substr(8);
            }
        }

        // clear buffer
        if(buf.size() > 0) {
            while(buf.size() < 8) {
                buf += "0";
            }
            byte temp = (byte)strtoul(buf.c_str(),NULL,2);
            out.write((char*)&temp,sizeof(temp));
        }
    }

    out.close();
    return 0;
}

/* Decompress the input files (regular Huffman decompression on each)
 * INPUT:  The prefix of the files to decompress
 * OUTPUT: The uncompressed file
 */
int decompress( string INFILE, string OUTFILE ) {
    // set up IO streams
    ifstream in(INFILE, ios::binary);
    ofstream out(OUTFILE);

    // decompress file
    while(true) {
        byte top = in.get();
        if(top == 255) {
            break;
        }
        int numChars;
        in.read((char*)&numChars,sizeof(numChars));
        if(top < 5) {
            char symbol = 'Z';
            switch(top) {
                case 0: symbol = 'A'; break;
                case 1: symbol = 'C'; break;
                case 2: symbol = 'G'; break;
                case 3: symbol = 'T'; break;
                case 4: symbol = 'N'; break;
                default: cerr << "ERROR: Unrecognized topology: " << top << endl << endl; exit(-1);
            }
            for(int i = 0; i < numChars; ++i) {
                out << symbol;
            }
        }
        else {
            Node* root = buildTree((int)top);
            Node* c = root;
            int printed = 0;
            while(printed < numChars) {
                byte buf = in.get();
                for(int i = 7; i >= 0; --i) {
                    int bit = (buf >> i) & 1;
                    if(bit == 0) {
                        if(c->c0 != NULL) {
                            c = c->c0;
                            if(c->c1 == NULL && c->c0 == NULL) {
                                out << c->symbol;
                                if(++printed == numChars) {
                                    break;
                                }
                                c = root;
                            }
                        }
                        else {
                            cerr << "ERROR: Invalid c0" << endl << endl; exit(-1);
                        }
                    }
                    else {
                        if(c->c1 != NULL) {
                            c = c->c1;
                            if(c->c1 == NULL && c->c0 == NULL) {
                                out << c->symbol;
                                if(++printed == numChars) {
                                    break;
                                }
                                c = root;
                            }
                        }
                        else {
                            cerr << "ERROR: Invalid c1" << endl << endl; exit(-1);
                        }
                    }
                }
            }
        }
    }
    return 0;
}

/* MAIN METHOD
 */
int main( int argc, char* argv[] ) {
    // parse arguments
    if(argc != 3) {
        cerr << "ERROR: Incorrect number of arguments" << endl;
        cerr << "See file header for usage information" << endl << endl;
        exit(-1);
    }
    const string COMMAND = argv[1];
    const string IN = argv[2];
    if(COMMAND != "compress" && COMMAND != "decompress") {
        cerr << "ERROR: First argument must be \"compress\" or \"decompress\"!" << endl;
        cerr << "See file header for usage information" << endl << endl;
        exit(-1);
    }

    // get tree topologies
    for(int i = 0; i < NUMTOPS; ++i) {
        TOPS[i] = getCode(i);
    }

    // run relevant function
    if(COMMAND == "compress") {
        return compress(IN,IN+".hsf");
    }
    else {
        return decompress(IN,IN.substr(0,IN.find(".hsf")));
    }
}
