/* AUTHOR: Niema Moshiri
 * DNA Huffman Compress algorithm
 *
 * USAGE:
 * -Huffman:
 * --Huffman Compress:   java DnaHuffman compress <in_file> <out_file>
 * --Huffman Decompress: java DnaHuffman decompress <in_file> <out_file>
 * 
 * -SplitMerge:
 * --Split Compress:     java DnaHuffman splitcompress <in_file> <out_prefix>
 * --Merge Decompress:   java DnaHuffman mergedecompress <in_prefix> <out_file>
 *
 * COMPRESSED FILE OUTPUT FORMAT:
 * -The first byte of the compressed file ("InfoByte") tells us the tree topology:
 * --If the first byte is 00000000, then the topology is balanced, so use A = 00, C = 01, G = 10, T = 11
 * --Otherwise, the first two bits say how many unique symbols appear: 01XXXXXX = 1 unique symbol, etc.
 * --The remaining 6 bits list the nucleotides in descending frequency order, using 00 = A, 01 = C, 10 = G, and 11 = T
 * --If there are three unique symbols (i.e., 11XXXXXX), we infer the lowest frequency symbol by listing the first 3 in the 6 remaining bits
 * -The next 8 bytes of the compressed file represent a long telling us how many total symbols are in the compressed file
 * -The remaining bytes are the compressed data
 *
 * If there is only 1 unique symbol, the resulting compressed file will only contain the first 9 bytes ("InfoByte" + "numChars")
 */
import java.io.*;
import java.nio.file.Paths;
import java.nio.file.Files;
import java.util.*;

public class DnaHuffman {
    /* MAIN METHOD */
    public static void main( String[] args ) {
        if(args.length != 3) {
            System.err.println("ERROR: Incorrect number of arguments");
            System.err.println("See file header for usage information");
            System.exit(-1);
        }
        final String IN = args[1];
        final String OUT = args[2];
        int c = -1;
        switch(args[0]) {
            case "compress": compress(IN,OUT); break;
            case "decompress": decompress(IN,OUT); break;
            case "splitcompress": splitcompress(IN,OUT); break;
            case "mergedecompress": mergedecompress(IN,OUT); break;
            default: System.err.println("ERROR: First argument must be \"compress\", \"decompress\", \"splitcompress\", or \"mergedecompress\"!"); System.exit(-1);
        }
    }
    
    /* SPLIT COMPRESS (only supports strings less than ~1 GB in size)*/
    public static void splitcompress(String INFILE, String OUTPREFIX) {
        // read number of characters
        int numChars = 0;
        try {
            InputStream in = new BufferedInputStream(new FileInputStream(new File(INFILE)));
            int b;
            while((b = in.read()) != -1) {
                ++numChars;
            }
            in.close();
        }
        catch(FileNotFoundException e) {
            System.err.println("ERROR: File \"" + INFILE + "\" not found!");
            System.exit(-1);
        }
        catch(IOException e) {
            System.err.println("ERROR: IOException while reading \"" + INFILE + "\"!");
            System.exit(-1);
        }
        
        // construct F[][], where F[n][i] is the number of occurrences of nucleotide n before and including position i of INFILE. 0 = A, 1 = C, 2 = G, 3 = T
        int[][] F = new int[4][numChars+1];
        try {
            InputStream in = new BufferedInputStream(new FileInputStream(new File(INFILE)));
            int b;
            int i = 1;
            while((b = in.read()) != -1) {
                F[0][i] = F[0][i-1];
                F[1][i] = F[1][i-1];
                F[2][i] = F[2][i-1];
                F[3][i] = F[3][i-1];
                switch((char)b) {
                    case 'A': ++F[0][i]; break;
                    case 'C': ++F[1][i]; break;
                    case 'G': ++F[2][i]; break;
                    case 'T': ++F[3][i]; break;
                    default: System.err.println("ERROR: Invalid character: '" + (char)b + "'!"); System.exit(-1);
                }
                ++i;
            }
            in.close();
        }
        catch(FileNotFoundException e) {
            System.err.println("ERROR: File \"" + INFILE + "\" not found!");
            System.exit(-1);
        }
        catch(IOException e) {
            System.err.println("ERROR: IOException while reading \"" + INFILE + "\"!");
            System.exit(-1);
        }
        
        // find optimal cuts (include 0 and |Genome| as "cuts")
        //int[] cuts = {0,100,200,300,numChars}; // THIS IS A DUMMY TEST!!!!! FIX THIS!!!!!!!!! ONCE I SOLVE THE ALGORITHM
        int[] cuts = getOptimalCuts(F);
        String numfiles = ""+(cuts.length-1);
        //F = null; // free up space now that I've found cuts
        
        // read file into memory
        try {
            String in = new String(Files.readAllBytes(Paths.get(INFILE)));
            for(int c = 1; c < cuts.length; ++c) {
                int start = cuts[c-1];
                int end = cuts[c];
                String s = in.substring(start,end);
                
                // build Huffman Tree
                HuffNode[] leaves = new HuffNode[4]; // 0 = A, 1 = C, 2 = G, 3 = T
                leaves[0] = new HuffNode('A');
                leaves[1] = new HuffNode('C');
                leaves[2] = new HuffNode('G');
                leaves[3] = new HuffNode('T');
                leaves[0].count = F[0][end] - F[0][start];
                leaves[1].count = F[1][end] - F[1][start];
                leaves[2].count = F[2][end] - F[2][start];
                leaves[3].count = F[3][end] - F[3][start];
                PriorityQueue<HuffNode> pq = new PriorityQueue<>(4);
                for(HuffNode l : leaves) {
                    pq.add(l);
                }
                while(pq.size() != 1) {
                    HuffNode n1 = pq.poll();
                    HuffNode n2 = pq.poll();
                    HuffNode n3 = new HuffNode(n1.symbol);
                    if(n2.symbol < n3.symbol) {
                        n3.symbol = n2.symbol;
                    }
                    n3.count = n1.count + n2.count;
                    n3.l = n1;
                    n3.r = n2;
                    n1.p = n3;
                    n2.p = n3;
                    pq.add(n3);
                }
                HuffNode root = pq.poll();
                
                // find topology (assuming all 4 nucleotides appear)
                int numUnique = 0;
                int topology = 1; // topology 1 is balanced, topology 2 is unbalanced
                for(HuffNode l : leaves) {
                    if(l.count != 0) {
                        ++numUnique;
                    }
                    if(l.code().length() == 3) {
                        topology = 2;
                    }
                }
                
                // set infoByte and "nucleotide to binary" mapping
                String infoByte = "";
                HashMap<Character,String> nucToBin = new HashMap<>();
                if(numUnique == 1) { // one unique symbol, so map it to 1
                    char nuc = 'Z';
                    for(HuffNode l : leaves) {
                        if(l.count != 0) {
                            nuc = l.symbol;
                            break;
                        }
                    }
                    nucToBin.put(nuc,"1");
                    infoByte = "01"; // 1 unique symbol
                    switch(nuc) {    // which symbol is the 1 unique symbol?
                        case 'A': infoByte += "00"; break;
                        case 'C': infoByte += "01"; break;
                        case 'G': infoByte += "10"; break;
                        case 'T': infoByte += "11"; break;
                        default: System.exit(-1);
                    }
                    infoByte += "0000"; // need it to be 8 bits
                }
                else if(numUnique == 2) { // two unique symbols, so map larger to 1 and smaller to 0
                    HuffNode n1 = null;
                    HuffNode n2 = null;
                    for(HuffNode l : leaves) {
                        if(n1 == null || l.count > n1.count) {
                            n2 = n1;
                            n1 = l;
                        }
                        else if(n2 == null || l.count > n2.count) {
                            n2 = l;
                        }
                    }
                    nucToBin.put(n1.symbol,"1");
                    nucToBin.put(n2.symbol,"0");
                    infoByte = "10";    // 2 unique symbols
                    switch(n1.symbol) { // which symbol is larger?
                        case 'A': infoByte += "00"; break;
                        case 'C': infoByte += "01"; break;
                        case 'G': infoByte += "10"; break;
                        case 'T': infoByte += "11"; break;
                        default: System.exit(-1);
                    }
                    switch(n2.symbol) { // which symbol is smaller?
                        case 'A': infoByte += "00"; break;
                        case 'C': infoByte += "01"; break;
                        case 'G': infoByte += "10"; break;
                        case 'T': infoByte += "11"; break;
                        default: System.exit(-1);
                    }
                    infoByte += "00";  // need it to be 8 bits
                }
                else if(numUnique == 3) { // three unique symbols, so map largest to 1, next largest to 01, and smallest to 00
                    HuffNode n1 = null;
                    HuffNode n2 = null;
                    HuffNode n3 = null;
                    for(HuffNode l : leaves) {
                        if(n1 == null || l.count > n1.count) {
                            n3 = n2;
                            n2 = n1;
                            n1 = l;
                        }
                        else if(n2 == null || l.count > n2.count) {
                            n3 = n2;
                            n2 = l;
                        }
                        else if(n3 == null || l.count > n3.count) {
                            n3 = l;
                        }
                    }
                    nucToBin.put(n1.symbol,"1");
                    nucToBin.put(n2.symbol,"01");
                    nucToBin.put(n3.symbol,"00");
                    infoByte = "11";    // 3 unique symbols
                    switch(n1.symbol) { // which symbol is largest?
                        case 'A': infoByte += "00"; break;
                        case 'C': infoByte += "01"; break;
                        case 'G': infoByte += "10"; break;
                        case 'T': infoByte += "11"; break;
                        default: System.exit(-1);
                    }
                    switch(n2.symbol) { // which symbol is 2nd largest?
                        case 'A': infoByte += "00"; break;
                        case 'C': infoByte += "01"; break;
                        case 'G': infoByte += "10"; break;
                        case 'T': infoByte += "11"; break;
                        default: System.exit(-1);
                    }
                    switch(n3.symbol) { // which symbol is smallest?
                        case 'A': infoByte += "00"; break;
                        case 'C': infoByte += "01"; break;
                        case 'G': infoByte += "10"; break;
                        case 'T': infoByte += "11"; break;
                        default: System.exit(-1);
                    }
                }
                else if(topology == 1) { // 4 unique symbols, balanced topology (A = 00, C = 01, G = 10, T = 11)
                    infoByte = "00000000";
                    nucToBin = new HashMap<>();
                    nucToBin.put('A',"00");
                    nucToBin.put('C',"01");
                    nucToBin.put('G',"10");
                    nucToBin.put('T',"11");
                }
                else { // 4 unique symbols, unbalanced topology
                    HuffNode n1 = null;
                    HuffNode n2 = null;
                    HuffNode n3 = null;
                    HuffNode n4 = null;
                    for(HuffNode l : leaves) {
                        if(n1 == null || l.count > n1.count) {
                            n4 = n3;
                            n3 = n2;
                            n2 = n1;
                            n1 = l;
                        }
                        else if(n2 == null || l.count > n2.count) {
                            n4 = n3;
                            n3 = n2;
                            n2 = l;
                        }
                        else if(n3 == null || l.count > n3.count) {
                            n4 = n3;
                            n3 = l;
                        }
                        else {
                            n4 = l;
                        }
                    }
                    nucToBin.put(n1.symbol,"1");
                    nucToBin.put(n2.symbol,"01");
                    nucToBin.put(n3.symbol,"001");
                    nucToBin.put(n4.symbol,"000");
                    infoByte = "00";    // 4 unique symbols
                    switch(n1.symbol) { // which symbol is largest?
                        case 'A': infoByte += "00"; break;
                        case 'C': infoByte += "01"; break;
                        case 'G': infoByte += "10"; break;
                        case 'T': infoByte += "11"; break;
                        default: System.exit(-1);
                    }
                    switch(n2.symbol) { // which symbol is 2nd largest?
                        case 'A': infoByte += "00"; break;
                        case 'C': infoByte += "01"; break;
                        case 'G': infoByte += "10"; break;
                        case 'T': infoByte += "11"; break;
                        default: System.exit(-1);
                    }
                    switch(n3.symbol) { // which symbol is 3rd largest? smallest can be inferred from these 3
                        case 'A': infoByte += "00"; break;
                        case 'C': infoByte += "01"; break;
                        case 'G': infoByte += "10"; break;
                        case 'T': infoByte += "11"; break;
                        default: System.exit(-1);
                    }
                }
                
                // encode file
                String filenum = ""+c;
                while(filenum.length() < numfiles.length()) {
                    filenum = "0" + filenum;
                }
                DataOutputStream out = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(new File(OUTPREFIX + filenum))));
                
                // write infoByte, then numChars
                out.writeByte(Integer.parseInt(infoByte, 2));
                out.writeLong((long)s.length());
                
                // if one unique symbol, we're done with just infobyte and numchars
                if(numUnique == 1) {
                    out.close();
                    continue;
                }
                
                String buf = "";
                for(int i = 0; i < s.length(); ++i) {
                    char b = s.charAt(i);
                    buf += nucToBin.get(b);
                    if(buf.length() >= 8) {
                        out.writeByte(Integer.parseInt(buf.substring(0,8),2));
                        buf = buf.substring(8);
                    }
                }
                if(buf.length() > 0) {
                    while(buf.length() < 8) {
                        buf += "0";
                    }
                    out.writeByte(Integer.parseInt(buf,2));
                }
                out.close();
            }
        }
        catch(Exception e) {
            e.printStackTrace();
            System.exit(-1);
        }
    }
    
    /* GET OPTIMAL CUTS
       1 UNIQUE CHAR:
        0:  (A)
        1:  (C)
        2:  (G)
        3:  (T)
       2 UNIQUE CHARS:
        4:  (A,C)
        5:  (A,G)
        6:  (A,T)
        7:  (C,G)
        8:  (C,T)
        9:  (G,T)
       3 UNIQUE CHARS:
        10: (A,(C,G))
        11: (A,(C,T))
        12: (A,(G,T))
        13: (C,(A,G))
        14: (C,(A,T))
        15: (C,(G,T))
        16: (G,(A,C))
        17: (G,(A,T))
        18: (G,(C,T))
        19: (T,(A,C))
        20: (T,(A,G))
        21: (T,(C,G))
       4 UNIQUE CHARS (BALANCED):
        22: ((A,C),(G,T))
       4 UNIQUE CASES (UNBALANCED):
        23: (A,(C,(G,T)))
        24: (A,(G,(C,T)))
        25: (A,(T,(C,G)))
        26: (C,(A,(G,T)))
        27: (C,(G,(A,T)))
        28: (C,(T,(A,G)))
        29: (G,(A,(C,T)))
        30: (G,(C,(A,T)))
        31: (G,(T,(A,C)))
        32: (T,(A,(C,G)))
        33: (T,(C,(A,G)))
        34: (T,(G,(A,C)))
        
        F[][], where F[n][i] is the number of occurrences of nucleotide n before and including position i of INFILE. 0 = A, 1 = C, 2 = G, 3 = T
    */
    public static int[] getOptimalCuts(int[][] F) {
        // initialize DP matrix + backtrack
        int length = F[0].length-1; // length of string
        long inf = 80*length; // didn't want to use Integer.MAX_VALUE to avoid overflow
        long[][] C = new long[length][35];
        int[][] backtrack = new int[length][35];
        if(F[0][1] != F[0][0]) {      // first char is 'A'
            for(int t = 0; t < 35; ++t) {
                if(t == 1 || t == 2 || t == 3 || t == 7 || t == 8 || t == 9 || t == 15 || t == 18 || t == 21) {
                    C[0][t] = inf;
                    backtrack[0][t] = -1;
                }
                else {
                    C[0][t] = 72L;
                    if(t == 4 || t == 5 || t == 6 || t == 10 || t == 11 || t == 12 || t == 23 || t == 24 || t == 25) {
                        ++C[0][t];
                    }
                    else if(t == 13 || t == 14 || t == 16 || t == 17 || t == 19 || t == 20 || t == 22 || t == 26 || t == 29 || t == 32) {
                        C[0][t] += 2;
                    }
                    else if(t != 0) {
                        C[0][t] += 3;
                    }
                    backtrack[0][t] = t;
                }
            }
        }
        else if(F[1][1] != F[1][0]) { // first char is 'C'
            for(int t = 0; t < 35; ++t) {
                if(t == 0 || t == 2 || t == 3 || t == 5 || t == 6 || t == 9 || t == 12 || t == 17 || t == 20) {
                    C[0][t] = inf;
                    backtrack[0][t] = -1;
                }
                else {
                    C[0][t] = 72L;
                    if(t == 4 || t == 7 || t == 8 || t == 13 || t == 14 || t == 15 || t == 26 || t == 27 || t == 28) {
                        ++C[0][t];
                    }
                    else if(t == 10 || t == 11 || t == 16 || t == 18 || t == 19 || t == 21 || t == 22 || t == 23 || t == 30 || t == 33 || t == 34) {
                        C[0][t] += 2;
                    }
                    else if(t != 1) {
                        C[0][t] += 3;
                    }
                    backtrack[0][t] = t;
                }
            }
        }
        else if(F[2][1] != F[2][0]) { // first char is 'G'
            for(int t = 0; t < 35; ++t) {
                if(t == 0 || t == 1 || t == 3 || t == 4 || t == 6 || t == 8 || t == 11 || t == 14 || t == 19) {
                    C[0][t] = inf;
                    backtrack[0][t] = -1;
                }
                else {
                    C[0][t] = 72L;
                    if(t == 5 || t == 7 || t == 9 || t == 16 || t == 17 || t == 18 || t == 29 || t == 30 || t == 31) {
                        ++C[0][t];
                    }
                    else if(t == 10 || t == 12 || t == 13 || t == 15 || t == 20 || t == 21 || t == 22 || t == 24 || t == 27 || t == 34) {
                        C[0][t] += 2;
                    }
                    else if(t != 2)  {
                        C[0][t] += 3;
                    }
                    backtrack[0][t] = t;
                }
            }
        }
        else {                          // first char is 'T'
            for(int t = 0; t < 35; ++t) {
                if(t == 0 || t == 1 || t == 2 || t == 4 || t == 5 || t == 7 || t == 10 || t == 13 || t == 16) {
                    C[0][t] = inf;
                    backtrack[0][t] = -1;
                }
                else {
                    C[0][t] = 72L;
                    if(t == 6 || t == 8 || t == 9 || t == 19 || t == 20 || t == 21 || t == 32 || t == 33 || t == 34) {
                        ++C[0][t];
                    }
                    else if(t == 11 || t == 12 || t == 14 || t == 15 || t == 17 || t == 18 || t == 22 || t == 25 || t == 28 || t == 31) {
                        C[0][t] += 2;
                    }
                    else if(t != 3)  {
                        C[0][t] += 3;
                    }
                    backtrack[0][t] = t;
                }
            }
        }
        
        // solve DP
        for(int i = 1; i < length; ++i) {
            for(int top = 0; top < 35; ++top) {
                long bestC = inf;
                int bestT = -1;
                
                // get next char of message
                char next = (char)0;
                if(F[0][i+1] != F[0][i]) {
                    next = 'A';
                }
                else if(F[1][i+1] != F[1][i]) {
                    next = 'C';
                }
                else if(F[2][i+1] != F[2][i]) {
                    next = 'G';
                }
                else {
                    next = 'T';
                }
                
                // Topology 0:  (A)
                if(top == 0) {
                    if(next != 'A') {
                        C[i][top] = inf;
                        backtrack[i][top] = -1;
                    }
                    else {
                        for(int t = 0; t < 35; ++t) {
                            long cost = C[i-1][t];
                            if(t != top) {
                                cost += 72L;
                            }
                            if(cost < bestC && backtrack[i-1][t] != -1) {
                                bestC = cost;
                                bestT = t;
                            }
                        }
                    }
                }
                
                // Topology 1:  (C)
                else if(top == 1) {
                    if(next != 'C') {
                        C[i][top] = inf;
                        backtrack[i][top] = -1;
                    }
                    else {
                        for(int t = 0; t < 35; ++t) {
                            long cost = C[i-1][t];
                            if(t != top) {
                                cost += 72L;
                            }
                            if(cost < bestC && backtrack[i-1][t] != -1) {
                                bestC = cost;
                                bestT = t;
                            }
                        }
                    }
                }
                
                // Topology 2:  (G)
                else if(top == 2) {
                    if(next != 'G') {
                        C[i][top] = inf;
                        backtrack[i][top] = -1;
                    }
                    else {
                        for(int t = 0; t < 35; ++t) {
                            long cost = C[i-1][t];
                            if(t != top) {
                                cost += 72L;
                            }
                            if(cost < bestC && backtrack[i-1][t] != -1) {
                                bestC = cost;
                                bestT = t;
                            }
                        }
                    }
                }
                
                // Topology 3:  (T)
                else if(top == 3) {
                    if(next != 'T') {
                        C[i][top] = inf;
                        backtrack[i][top] = -1;
                    }
                    else {
                        for(int t = 0; t < 35; ++t) {
                            long cost = C[i-1][t];
                            if(t != top) {
                                cost += 72L;
                            }
                            if(cost < bestC && backtrack[i-1][t] != -1) {
                                bestC = cost;
                                bestT = t;
                            }
                        }
                    }
                }
                
                // Topology 4:  (A,C)
                else if(top == 4) {
                    if(next == 'G' || next == 'T') {
                        C[i][top] = inf;
                        backtrack[i][top] = -1;
                    }
                    else {
                        for(int t = 0; t < 35; ++t) {
                            long cost = C[i-1][t] + 1;
                            if(t != top) {
                                cost += 72L;
                            }
                            if(cost < bestC && backtrack[i-1][t] != -1) {
                                bestC = cost;
                                bestT = t;
                            }
                        }
                    }
                }
                
                // Topology 5:  (A,G)
                else if(top == 5) {
                    if(next == 'C' || next == 'T') {
                        C[i][top] = inf;
                        backtrack[i][top] = -1;
                    }
                    else {
                        for(int t = 0; t < 35; ++t) {
                            long cost = C[i-1][t] + 1;
                            if(t != top) {
                                cost += 72L;
                            }
                            if(cost < bestC && backtrack[i-1][t] != -1) {
                                bestC = cost;
                                bestT = t;
                            }
                        }
                    }
                }
                
                // Topology 6:  (A,T)
                else if(top == 6) {
                    if(next == 'C' || next == 'G') {
                        C[i][top] = inf;
                        backtrack[i][top] = -1;
                    }
                    else {
                        for(int t = 0; t < 35; ++t) {
                            long cost = C[i-1][t] + 1;
                            if(t != top) {
                                cost += 72L;
                            }
                            if(cost < bestC && backtrack[i-1][t] != -1) {
                                bestC = cost;
                                bestT = t;
                            }
                        }
                    }
                }
                
                // Topology 7:  (C,G)
                else if(top == 7) {
                    if(next == 'A' || next == 'T') {
                        C[i][top] = inf;
                        backtrack[i][top] = -1;
                    }
                    else {
                        for(int t = 0; t < 35; ++t) {
                            long cost = C[i-1][t] + 1;
                            if(t != top) {
                                cost += 72L;
                            }
                            if(cost < bestC && backtrack[i-1][t] != -1) {
                                bestC = cost;
                                bestT = t;
                            }
                        }
                    }
                }
                
                // Topology 8:  (C,T)
                else if(top == 8) {
                    if(next == 'A' || next == 'G') {
                        C[i][top] = inf;
                        backtrack[i][top] = -1;
                    }
                    else {
                        for(int t = 0; t < 35; ++t) {
                            long cost = C[i-1][t] + 1;
                            if(t != top) {
                                cost += 72L;
                            }
                            if(cost < bestC && backtrack[i-1][t] != -1) {
                                bestC = cost;
                                bestT = t;
                            }
                        }
                    }
                }
                
                // Topology 9:  (G,T)
                else if(top == 9) {
                    if(next == 'A' || next == 'C') {
                        C[i][top] = inf;
                        backtrack[i][top] = -1;
                    }
                    else {
                        for(int t = 0; t < 35; ++t) {
                            long cost = C[i-1][t] + 1;
                            if(t != top) {
                                cost += 72L;
                            }
                            if(cost < bestC && backtrack[i-1][t] != -1) {
                                bestC = cost;
                                bestT = t;
                            }
                        }
                    }
                }
                
                // Topology 10: (A,(C,G))
                else if(top == 10) {
                    if(next == 'T') {
                        C[i][top] = inf;
                        backtrack[i][top] = -1;
                    }
                    else {
                        for(int t = 0; t < 35; ++t) {
                            long cost = C[i-1][t];
                            switch(next) {
                                case 'A': cost += 1; break;
                                case 'C': cost += 2; break;
                                case 'G': cost += 2; break;
                                default: System.err.println("SOMETHING WENT WRONG"); System.exit(-1);
                            }
                            if(t != top) {
                                cost += 72L;
                            }
                            if(cost < bestC && backtrack[i-1][t] != -1) {
                                bestC = cost;
                                bestT = t;
                            }
                        }
                    }
                }
                
                // Topology 11: (A,(C,T))
                else if(top == 11) {
                    if(next == 'G') {
                        C[i][top] = inf;
                        backtrack[i][top] = -1;
                    }
                    else {
                        for(int t = 0; t < 35; ++t) {
                            long cost = C[i-1][t];
                            switch(next) {
                                case 'A': cost += 1; break;
                                case 'C': cost += 2; break;
                                case 'T': cost += 2; break;
                                default: System.err.println("SOMETHING WENT WRONG"); System.exit(-1);
                            }
                            if(t != top) {
                                cost += 72L;
                            }
                            if(cost < bestC && backtrack[i-1][t] != -1) {
                                bestC = cost;
                                bestT = t;
                            }
                        }
                    }
                }
                
                // Topology 12: (A,(G,T))
                else if(top == 12) {
                    if(next == 'C') {
                        C[i][top] = inf;
                        backtrack[i][top] = -1;
                    }
                    else {
                        for(int t = 0; t < 35; ++t) {
                            long cost = C[i-1][t];
                            switch(next) {
                                case 'A': cost += 1; break;
                                case 'G': cost += 2; break;
                                case 'T': cost += 2; break;
                                default: System.err.println("SOMETHING WENT WRONG"); System.exit(-1);
                            }
                            if(t != top) {
                                cost += 72L;
                            }
                            if(cost < bestC && backtrack[i-1][t] != -1) {
                                bestC = cost;
                                bestT = t;
                            }
                        }
                    }
                }
                
                // Topology 13: (C,(A,G))
                else if(top == 13) {
                    if(next == 'T') {
                        C[i][top] = inf;
                        backtrack[i][top] = -1;
                    }
                    else {
                        for(int t = 0; t < 35; ++t) {
                            long cost = C[i-1][t];
                            switch(next) {
                                case 'C': cost += 1; break;
                                case 'A': cost += 2; break;
                                case 'G': cost += 2; break;
                                default: System.err.println("SOMETHING WENT WRONG"); System.exit(-1);
                            }
                            if(t != top) {
                                cost += 72L;
                            }
                            if(cost < bestC && backtrack[i-1][t] != -1) {
                                bestC = cost;
                                bestT = t;
                            }
                        }
                    }
                }
                
                // Topology 14: (C,(A,T))
                else if(top == 14) {
                    if(next == 'G') {
                        C[i][top] = inf;
                        backtrack[i][top] = -1;
                    }
                    else {
                        for(int t = 0; t < 35; ++t) {
                            long cost = C[i-1][t];
                            switch(next) {
                                case 'C': cost += 1; break;
                                case 'A': cost += 2; break;
                                case 'T': cost += 2; break;
                                default: System.err.println("SOMETHING WENT WRONG"); System.exit(-1);
                            }
                            if(t != top) {
                                cost += 72L;
                            }
                            if(cost < bestC && backtrack[i-1][t] != -1) {
                                bestC = cost;
                                bestT = t;
                            }
                        }
                    }
                }
                
                // Topology 15: (C,(G,T))
                else if(top == 15) {
                    if(next == 'A') {
                        C[i][top] = inf;
                        backtrack[i][top] = -1;
                    }
                    else {
                        for(int t = 0; t < 35; ++t) {
                            long cost = C[i-1][t];
                            switch(next) {
                                case 'C': cost += 1; break;
                                case 'G': cost += 2; break;
                                case 'T': cost += 2; break;
                                default: System.err.println("SOMETHING WENT WRONG"); System.exit(-1);
                            }
                            if(t != top) {
                                cost += 72L;
                            }
                            if(cost < bestC && backtrack[i-1][t] != -1) {
                                bestC = cost;
                                bestT = t;
                            }
                        }
                    }
                }
                
                // Topology 16: (G,(A,C))
                else if(top == 16) {
                    if(next == 'T') {
                        C[i][top] = inf;
                        backtrack[i][top] = -1;
                    }
                    else {
                        for(int t = 0; t < 35; ++t) {
                            long cost = C[i-1][t];
                            switch(next) {
                                case 'G': cost += 1; break;
                                case 'A': cost += 2; break;
                                case 'C': cost += 2; break;
                                default: System.err.println("SOMETHING WENT WRONG"); System.exit(-1);
                            }
                            if(t != top) {
                                cost += 72L;
                            }
                            if(cost < bestC && backtrack[i-1][t] != -1) {
                                bestC = cost;
                                bestT = t;
                            }
                        }
                    }
                }
                
                // Topology 17: (G,(A,T))
                else if(top == 17) {
                    if(next == 'C') {
                        C[i][top] = inf;
                        backtrack[i][top] = -1;
                    }
                    else {
                        for(int t = 0; t < 35; ++t) {
                            long cost = C[i-1][t];
                            switch(next) {
                                case 'G': cost += 1; break;
                                case 'A': cost += 2; break;
                                case 'T': cost += 2; break;
                                default: System.err.println("SOMETHING WENT WRONG"); System.exit(-1);
                            }
                            if(t != top) {
                                cost += 72L;
                            }
                            if(cost < bestC && backtrack[i-1][t] != -1) {
                                bestC = cost;
                                bestT = t;
                            }
                        }
                    }
                }
                
                // Topology 18: (G,(C,T))
                else if(top == 18) {
                    if(next == 'A') {
                        C[i][top] = inf;
                        backtrack[i][top] = -1;
                    }
                    else {
                        for(int t = 0; t < 35; ++t) {
                            long cost = C[i-1][t];
                            switch(next) {
                                case 'G': cost += 1; break;
                                case 'C': cost += 2; break;
                                case 'T': cost += 2; break;
                                default: System.err.println("SOMETHING WENT WRONG"); System.exit(-1);
                            }
                            if(t != top) {
                                cost += 72L;
                            }
                            if(cost < bestC && backtrack[i-1][t] != -1) {
                                bestC = cost;
                                bestT = t;
                            }
                        }
                    }
                }
                
                // Topology 19: (T,(A,C))
                else if(top == 19) {
                    if(next == 'G') {
                        C[i][top] = inf;
                        backtrack[i][top] = -1;
                    }
                    else {
                        for(int t = 0; t < 35; ++t) {
                            long cost = C[i-1][t];
                            switch(next) {
                                case 'T': cost += 1; break;
                                case 'A': cost += 2; break;
                                case 'C': cost += 2; break;
                                default: System.err.println("SOMETHING WENT WRONG"); System.exit(-1);
                            }
                            if(t != top) {
                                cost += 72L;
                            }
                            if(cost < bestC && backtrack[i-1][t] != -1) {
                                bestC = cost;
                                bestT = t;
                            }
                        }
                    }
                }
                
                // Topology 20: (T,(A,G))
                else if(top == 20) {
                    if(next == 'C') {
                        C[i][top] = inf;
                        backtrack[i][top] = -1;
                    }
                    else {
                        for(int t = 0; t < 35; ++t) {
                            long cost = C[i-1][t];
                            switch(next) {
                                case 'T': cost += 1; break;
                                case 'A': cost += 2; break;
                                case 'G': cost += 2; break;
                                default: System.err.println("SOMETHING WENT WRONG"); System.exit(-1);
                            }
                            if(t != top) {
                                cost += 72L;
                            }
                            if(cost < bestC && backtrack[i-1][t] != -1) {
                                bestC = cost;
                                bestT = t;
                            }
                        }
                    }
                }
                
                // Topology 21: (T,(C,G))
                else if(top == 21) {
                    if(next == 'A') {
                        C[i][top] = inf;
                        backtrack[i][top] = -1;
                    }
                    else {
                        for(int t = 0; t < 35; ++t) {
                            long cost = C[i-1][t];
                            switch(next) {
                                case 'T': cost += 1; break;
                                case 'C': cost += 2; break;
                                case 'G': cost += 2; break;
                                default: System.err.println("SOMETHING WENT WRONG"); System.exit(-1);
                            }
                            if(t != top) {
                                cost += 72L;
                            }
                            if(cost < bestC && backtrack[i-1][t] != -1) {
                                bestC = cost;
                                bestT = t;
                            }
                        }
                    }
                }
                
                // Topology 22: ((A,C),(G,T))
                else if(top == 22) {
                    for(int t = 0; t < 35; ++t) {
                        long cost = C[i-1][t] + 2;
                        if(t != top) {
                            cost += 72L;
                        }
                        if(cost < bestC && backtrack[i-1][t] != -1) {
                            bestC = cost;
                            bestT = t;
                        }
                    }
                }
                
                // Topology 23: (A,(C,(G,T)))
                else if(top == 23) {
                    for(int t = 0; t < 35; ++t) {
                        long cost = C[i-1][t];
                        switch(next) {
                            case 'A': cost += 1; break;
                            case 'C': cost += 2; break;
                            case 'G': cost += 3; break;
                            case 'T': cost += 3; break;
                            default: System.err.println("SOMETHING WENT WRONG"); System.exit(-1);
                        }
                        if(t != top) {
                            cost += 72L;
                        }
                        if(cost < bestC && backtrack[i-1][t] != -1) {
                            bestC = cost;
                            bestT = t;
                        }
                    }
                }
                
                // Topology 24: (A,(G,(C,T)))
                else if(top == 24) {
                    for(int t = 0; t < 35; ++t) {
                        long cost = C[i-1][t];
                        switch(next) {
                            case 'A': cost += 1; break;
                            case 'G': cost += 2; break;
                            case 'C': cost += 3; break;
                            case 'T': cost += 3; break;
                            default: System.err.println("SOMETHING WENT WRONG"); System.exit(-1);
                        }
                        if(t != top) {
                            cost += 72L;
                        }
                        if(cost < bestC && backtrack[i-1][t] != -1) {
                            bestC = cost;
                            bestT = t;
                        }
                    }
                }
                
                // Topology 25: (A,(T,(C,G)))
                else if(top == 25) {
                    for(int t = 0; t < 35; ++t) {
                        long cost = C[i-1][t];
                        switch(next) {
                            case 'A': cost += 1; break;
                            case 'T': cost += 2; break;
                            case 'C': cost += 3; break;
                            case 'G': cost += 3; break;
                            default: System.err.println("SOMETHING WENT WRONG"); System.exit(-1);
                        }
                        if(t != top) {
                            cost += 72L;
                        }
                        if(cost < bestC && backtrack[i-1][t] != -1) {
                            bestC = cost;
                            bestT = t;
                        }
                    }
                }
                
                // Topology 26: (C,(A,(G,T)))
                else if(top == 26) {
                    for(int t = 0; t < 35; ++t) {
                        long cost = C[i-1][t];
                        switch(next) {
                            case 'C': cost += 1; break;
                            case 'A': cost += 2; break;
                            case 'G': cost += 3; break;
                            case 'T': cost += 3; break;
                            default: System.err.println("SOMETHING WENT WRONG"); System.exit(-1);
                        }
                        if(t != top) {
                            cost += 72L;
                        }
                        if(cost < bestC && backtrack[i-1][t] != -1) {
                            bestC = cost;
                            bestT = t;
                        }
                    }
                }
                
                // Topology 27: (C,(G,(A,T)))
                else if(top == 27) {
                    for(int t = 0; t < 35; ++t) {
                        long cost = C[i-1][t];
                        switch(next) {
                            case 'C': cost += 1; break;
                            case 'G': cost += 2; break;
                            case 'A': cost += 3; break;
                            case 'T': cost += 3; break;
                            default: System.err.println("SOMETHING WENT WRONG"); System.exit(-1);
                        }
                        if(t != top) {
                            cost += 72L;
                        }
                        if(cost < bestC && backtrack[i-1][t] != -1) {
                            bestC = cost;
                            bestT = t;
                        }
                    }
                }
                
                // Topology 28: (C,(T,(A,G)))
                else if(top == 28) {
                    for(int t = 0; t < 35; ++t) {
                        long cost = C[i-1][t];
                        switch(next) {
                            case 'C': cost += 1; break;
                            case 'T': cost += 2; break;
                            case 'A': cost += 3; break;
                            case 'G': cost += 3; break;
                            default: System.err.println("SOMETHING WENT WRONG"); System.exit(-1);
                        }
                        if(t != top) {
                            cost += 72L;
                        }
                        if(cost < bestC && backtrack[i-1][t] != -1) {
                            bestC = cost;
                            bestT = t;
                        }
                    }
                }
                
                // Topology 29: (G,(A,(C,T)))
                else if(top == 29) {
                    for(int t = 0; t < 35; ++t) {
                        long cost = C[i-1][t];
                        switch(next) {
                            case 'G': cost += 1; break;
                            case 'A': cost += 2; break;
                            case 'C': cost += 3; break;
                            case 'T': cost += 3; break;
                            default: System.err.println("SOMETHING WENT WRONG"); System.exit(-1);
                        }
                        if(t != top) {
                            cost += 72L;
                        }
                        if(cost < bestC && backtrack[i-1][t] != -1) {
                            bestC = cost;
                            bestT = t;
                        }
                    }
                }
                
                // Topology 30: (G,(C,(A,T)))
                else if(top == 30) {
                    for(int t = 0; t < 35; ++t) {
                        long cost = C[i-1][t];
                        switch(next) {
                            case 'G': cost += 1; break;
                            case 'C': cost += 2; break;
                            case 'A': cost += 3; break;
                            case 'T': cost += 3; break;
                            default: System.err.println("SOMETHING WENT WRONG"); System.exit(-1);
                        }
                        if(t != top) {
                            cost += 72L;
                        }
                        if(cost < bestC && backtrack[i-1][t] != -1) {
                            bestC = cost;
                            bestT = t;
                        }
                    }
                }
                
                // Topology 31: (G,(T,(A,C)))
                else if(top == 31) {
                    for(int t = 0; t < 35; ++t) {
                        long cost = C[i-1][t];
                        switch(next) {
                            case 'G': cost += 1; break;
                            case 'T': cost += 2; break;
                            case 'A': cost += 3; break;
                            case 'C': cost += 3; break;
                            default: System.err.println("SOMETHING WENT WRONG"); System.exit(-1);
                        }
                        if(t != top) {
                            cost += 72L;
                        }
                        if(cost < bestC && backtrack[i-1][t] != -1) {
                            bestC = cost;
                            bestT = t;
                        }
                    }
                }
                
                // Topology 32: (T,(A,(C,G)))
                else if(top == 32) {
                    for(int t = 0; t < 35; ++t) {
                        long cost = C[i-1][t];
                        switch(next) {
                            case 'T': cost += 1; break;
                            case 'A': cost += 2; break;
                            case 'C': cost += 3; break;
                            case 'G': cost += 3; break;
                            default: System.err.println("SOMETHING WENT WRONG"); System.exit(-1);
                        }
                        if(t != top) {
                            cost += 72L;
                        }
                        if(cost < bestC && backtrack[i-1][t] != -1) {
                            bestC = cost;
                            bestT = t;
                        }
                    }
                }
                
                // Topology 33: (T,(C,(A,G)))
                else if(top == 33) {
                    for(int t = 0; t < 35; ++t) {
                        long cost = C[i-1][t];
                        switch(next) {
                            case 'T': cost += 1; break;
                            case 'C': cost += 2; break;
                            case 'A': cost += 3; break;
                            case 'G': cost += 3; break;
                            default: System.err.println("SOMETHING WENT WRONG"); System.exit(-1);
                        }
                        if(t != top) {
                            cost += 72L;
                        }
                        if(cost < bestC && backtrack[i-1][t] != -1) {
                            bestC = cost;
                            bestT = t;
                        }
                    }
                }
                
                // Topology 34: (T,(G,(A,C)))
                else if(top == 34) {
                    for(int t = 0; t < 35; ++t) {
                        long cost = C[i-1][t];
                        switch(next) {
                            case 'T': cost += 1; break;
                            case 'G': cost += 2; break;
                            case 'A': cost += 3; break;
                            case 'C': cost += 3; break;
                            default: System.err.println("SOMETHING WENT WRONG"); System.exit(-1);
                        }
                        if(t != top) {
                            cost += 72L;
                        }
                        if(cost < bestC && backtrack[i-1][t] != -1) {
                            bestC = cost;
                            bestT = t;
                        }
                    }
                }
                
                // set values
                C[i][top] = bestC;
                backtrack[i][top] = bestT;
            }
        }
        
        // NIEMA TEST
        /*for(int row = 0; row < C.length; ++row) {
            for(int col = 0; col < C[0].length; ++col) {
                System.out.print(backtrack[row][col] + ",");
            }
            System.out.println();
        }*/
        
        // list of positions at which to make cuts
        ArrayList<Integer> cuts = new ArrayList<>();
        cuts.add(0);
        cuts.add(length);
        int top = 0;
        for(int t = 1; t < 35; ++t) {
            if(C[length-1][t] < C[length-1][top]) {
                top = t;
            }
        }
        for(int i = length-1; i > 0; --i) {
            if(top != backtrack[i][top]) {
                cuts.add(i);
            }
            top = backtrack[i][top];
        }
        int[] out = new int[cuts.size()];
        for(int i = 0; i < out.length; ++i) {
            out[i] = cuts.get(i);
        }
        Arrays.sort(out);
        return out;
    }
    
    /* MERGE DECOMPRESS */
    public static void mergedecompress(String INPREFIX, String OUTFILE) {
        try {
            // set up files
            DataOutputStream out = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(new File(OUTFILE))));
            File dir = new File(".");
            File[] foundFiles = dir.listFiles(new FilenameFilter() {
                public boolean accept(File dir, String name) {
                    return name.startsWith(INPREFIX);
                }
            });
            
            // decompress each file
            for(File f : foundFiles) {
                DataInputStream in = new DataInputStream(new BufferedInputStream(new FileInputStream(f)));
                byte infoByte = in.readByte();
                int numUnique = (int)(infoByte >>> 6) & 3;
                
                // if 1 unique symbol, just read length and print that symbol length times
                if(numUnique == 1) {
                    char symbol = 'Z';
                    switch((int)(infoByte >>> 4) & 3) { // the symbol
                        case 0: symbol = 'A'; break;
                        case 1: symbol = 'C'; break;
                        case 2: symbol = 'G'; break;
                        case 3: symbol = 'T'; break;
                        default: System.exit(-1);
                    }
                    long numChars = in.readLong();
                    for(long i = 0; i < numChars; ++i) {
                        out.writeByte((byte)symbol);
                    }
                    in.close();
                    continue;
                }
                
                // build Huffman Tree
                HuffNode root = new HuffNode('\0');
                if((int)infoByte == 0) { // all 4 symbols, topology 1, so balanced
                    root.l = new HuffNode('\0');
                    root.r = new HuffNode('\0');
                    root.r.r = new HuffNode('A'); // 00 = A
                    root.r.l = new HuffNode('C'); // 01 = C
                    root.l.r = new HuffNode('G'); // 10 = G
                    root.l.l = new HuffNode('T'); // 11 = T
                }
                else if(numUnique == 2) { // 2 unique symbols
                    switch((int)(infoByte >>> 4) & 3) { // the larger symbol
                        case 0: root.l = new HuffNode('A'); break;
                        case 1: root.l = new HuffNode('C'); break;
                        case 2: root.l = new HuffNode('G'); break;
                        case 3: root.l = new HuffNode('T'); break;
                        default: System.exit(-1);
                    }
                    switch((int)(infoByte >>> 2) & 3) { // the smaller symbol
                        case 0: root.r = new HuffNode('A'); break;
                        case 1: root.r = new HuffNode('C'); break;
                        case 2: root.r = new HuffNode('G'); break;
                        case 3: root.r = new HuffNode('T'); break;
                        default: System.exit(-1);
                    }
                }
                else if(numUnique == 3) { // 3 unique symbols
                    root.r = new HuffNode('\0');
                    switch((int)(infoByte >>> 4) & 3) { // the largest symbol
                        case 0: root.l = new HuffNode('A'); break;
                        case 1: root.l = new HuffNode('C'); break;
                        case 2: root.l = new HuffNode('G'); break;
                        case 3: root.l = new HuffNode('T'); break;
                        default: System.exit(-1);
                    }
                    switch((int)(infoByte >>> 2) & 3) { // the middle symbol
                        case 0: root.r.l = new HuffNode('A'); break;
                        case 1: root.r.l = new HuffNode('C'); break;
                        case 2: root.r.l = new HuffNode('G'); break;
                        case 3: root.r.l = new HuffNode('T'); break;
                        default: System.exit(-1);
                    }
                    switch((int)infoByte & 3) { // the smallest symbol
                        case 0: root.r.r = new HuffNode('A'); break;
                        case 1: root.r.r = new HuffNode('C'); break;
                        case 2: root.r.r = new HuffNode('G'); break;
                        case 3: root.r.r = new HuffNode('T'); break;
                        default: System.exit(-1);
                    }
                }
                else { // 4 unique symbols, unbalanced topology
                    char l1 = 'Z';
                    char l2 = 'Z';
                    char l3 = 'Z';
                    boolean[] used = new boolean[4];
                    switch((int)(infoByte >>> 4) & 3) { // extract l1
                        case 0: l1 = 'A'; used[0] = true; break;
                        case 1: l1 = 'C'; used[1] = true; break;
                        case 2: l1 = 'G'; used[2] = true; break;
                        case 3: l1 = 'T'; used[3] = true; break;
                        default: System.exit(-1);
                    }
                    switch((int)(infoByte >>> 2) & 3) { // extract l2
                        case 0: l2 = 'A'; used[0] = true; break;
                        case 1: l2 = 'C'; used[1] = true; break;
                        case 2: l2 = 'G'; used[2] = true; break;
                        case 3: l2 = 'T'; used[3] = true; break;
                        default: System.exit(-1);
                    }
                    switch((int)infoByte & 3) { // extract l3
                        case 0: l3 = 'A'; used[0] = true; break;
                        case 1: l3 = 'C'; used[1] = true; break;
                        case 2: l3 = 'G'; used[2] = true; break;
                        case 3: l3 = 'T'; used[3] = true; break;
                        default: System.exit(-1);
                    }
                    root.l = new HuffNode(l1); // length-1 symbol
                    root.r = new HuffNode('\0');
                    root.r.l = new HuffNode(l2); // length-2 symbol
                    root.r.r = new HuffNode('\0');
                    root.r.r.l = new HuffNode(l3); // length-3 symbol
                    for(int i = 0; i < 4; ++i) {
                        if(!used[i]) {
                            switch(i) {
                                case 0: root.r.r.r = new HuffNode('A'); break;
                                case 1: root.r.r.r = new HuffNode('C'); break;
                                case 2: root.r.r.r = new HuffNode('G'); break;
                                case 3: root.r.r.r = new HuffNode('T'); break;
                            }
                        }
                    }
                }
                
                // read rest of file and decode
                long numChars = in.readLong();
                long printed = 0L;
                HuffNode c = root;
                boolean done = false;
                while(!done) {
                    byte buf = in.readByte();
                    for(int i = 7; i >= 0; --i) {
                        int bit = (buf >>> i) & 1;
                        if(bit == 0) {
                            if(c.r != null) {
                                c = c.r;
                                if(c.l == null && c.r == null) {
                                    out.writeByte((byte)(c.symbol));
                                    if(++printed == numChars) {
                                        in.close();
                                        done = true;
                                        break;
                                    }
                                    c = root;
                                }
                            }
                            else {
                                c = root.r;
                            }
                        }
                        else {
                            if(c.l != null) {
                                c = c.l;
                                if(c.l == null && c.r == null) {
                                    out.writeByte((byte)(c.symbol));
                                    if(++printed == numChars) {
                                        in.close();
                                        done = true;
                                        break;
                                    }
                                    c = root;
                                }
                            }
                            else {
                                c = root.l;
                            }
                        }
                    }
                }
            }
            out.close();
            System.exit(0);
        }
        catch(Exception e) {
            e.printStackTrace();
            System.exit(-1);
        }
    }
    
    /* HUFFMAN COMPRESS */
    public static void compress(String INFILE, String OUTFILE) {
        // read file
        long numChars = 0L;
        HuffNode[] leaves = new HuffNode[4]; // 0 = A, 1 = C, 2 = G, 3 = T
        leaves[0] = new HuffNode('A');
        leaves[1] = new HuffNode('C');
        leaves[2] = new HuffNode('G');
        leaves[3] = new HuffNode('T');
        try {
            InputStream in = new BufferedInputStream(new FileInputStream(new File(INFILE)));
            int b;
            while((b = in.read()) != -1) {
                switch((char)b) {
                    case 'A': ++leaves[0].count; break;
                    case 'C': ++leaves[1].count; break;
                    case 'G': ++leaves[2].count; break;
                    case 'T': ++leaves[3].count; break;
                    default: System.err.println("ERROR: Invalid character: '" + (char)b + "'!"); System.exit(-1);
                }
                ++numChars;
            }
            in.close();
        }
        catch(FileNotFoundException e) {
            System.err.println("ERROR: File \"" + INFILE + "\" not found!");
            System.exit(-1);
        }
        catch(IOException e) {
            System.err.println("ERROR: IOException while reading \"" + INFILE + "\"!");
            System.exit(-1);
        }
        
        // build Huffman Tree
        PriorityQueue<HuffNode> pq = new PriorityQueue<>(4);
        for(HuffNode l : leaves) {
            pq.add(l);
        }
        while(pq.size() != 1) {
            HuffNode n1 = pq.poll();
            HuffNode n2 = pq.poll();
            HuffNode n3 = new HuffNode(n1.symbol);
            if(n2.symbol < n3.symbol) {
                n3.symbol = n2.symbol;
            }
            n3.count = n1.count + n2.count;
            n3.l = n1;
            n3.r = n2;
            n1.p = n3;
            n2.p = n3;
            pq.add(n3);
        }
        HuffNode root = pq.poll();
        
        // find topology (assuming all 4 nucleotides appear)
        int numUnique = 0;
        int topology = 1; // topology 1 is balanced, topology 2 is unbalanced
        for(HuffNode l : leaves) {
            if(l.count != 0) {
                ++numUnique;
            }
            if(l.code().length() == 3) {
                topology = 2;
            }
        }
        
        // set infoByte and "nucleotide to binary" mapping
        String infoByte = "";
        HashMap<Character,String> nucToBin = new HashMap<>();
        if(numUnique == 1) { // one unique symbol, so map it to 1
            char nuc = 'Z';
            for(HuffNode l : leaves) {
                if(l.count != 0) {
                    nuc = l.symbol;
                    break;
                }
            }
            nucToBin.put(nuc,"1");
            infoByte = "01"; // 1 unique symbol
            switch(nuc) {    // which symbol is the 1 unique symbol?
                case 'A': infoByte += "00"; break;
                case 'C': infoByte += "01"; break;
                case 'G': infoByte += "10"; break;
                case 'T': infoByte += "11"; break;
                default: System.exit(-1);
            }
            infoByte += "0000"; // need it to be 8 bits
        }
        else if(numUnique == 2) { // two unique symbols, so map larger to 1 and smaller to 0
            HuffNode n1 = null;
            HuffNode n2 = null;
            for(HuffNode l : leaves) {
                if(n1 == null || l.count > n1.count) {
                    n2 = n1;
                    n1 = l;
                }
                else if(n2 == null || l.count > n2.count) {
                    n2 = l;
                }
            }
            nucToBin.put(n1.symbol,"1");
            nucToBin.put(n2.symbol,"0");
            infoByte = "10";    // 2 unique symbols
            switch(n1.symbol) { // which symbol is larger?
                case 'A': infoByte += "00"; break;
                case 'C': infoByte += "01"; break;
                case 'G': infoByte += "10"; break;
                case 'T': infoByte += "11"; break;
                default: System.exit(-1);
            }
            switch(n2.symbol) { // which symbol is smaller?
                case 'A': infoByte += "00"; break;
                case 'C': infoByte += "01"; break;
                case 'G': infoByte += "10"; break;
                case 'T': infoByte += "11"; break;
                default: System.exit(-1);
            }
            infoByte += "00";  // need it to be 8 bits
        }
        else if(numUnique == 3) { // three unique symbols, so map largest to 1, next largest to 01, and smallest to 00
            HuffNode n1 = null;
            HuffNode n2 = null;
            HuffNode n3 = null;
            for(HuffNode l : leaves) {
                if(n1 == null || l.count > n1.count) {
                    n3 = n2;
                    n2 = n1;
                    n1 = l;
                }
                else if(n2 == null || l.count > n2.count) {
                    n3 = n2;
                    n2 = l;
                }
                else if(n3 == null || l.count > n3.count) {
                    n3 = l;
                }
            }
            nucToBin.put(n1.symbol,"1");
            nucToBin.put(n2.symbol,"01");
            nucToBin.put(n3.symbol,"00");
            infoByte = "11";    // 3 unique symbols
            switch(n1.symbol) { // which symbol is largest?
                case 'A': infoByte += "00"; break;
                case 'C': infoByte += "01"; break;
                case 'G': infoByte += "10"; break;
                case 'T': infoByte += "11"; break;
                default: System.exit(-1);
            }
            switch(n2.symbol) { // which symbol is 2nd largest?
                case 'A': infoByte += "00"; break;
                case 'C': infoByte += "01"; break;
                case 'G': infoByte += "10"; break;
                case 'T': infoByte += "11"; break;
                default: System.exit(-1);
            }
            switch(n3.symbol) { // which symbol is smallest?
                case 'A': infoByte += "00"; break;
                case 'C': infoByte += "01"; break;
                case 'G': infoByte += "10"; break;
                case 'T': infoByte += "11"; break;
                default: System.exit(-1);
            }
        }
        else if(topology == 1) { // 4 unique symbols, balanced topology (A = 00, C = 01, G = 10, T = 11)
            infoByte = "00000000";
            nucToBin = new HashMap<>();
            nucToBin.put('A',"00");
            nucToBin.put('C',"01");
            nucToBin.put('G',"10");
            nucToBin.put('T',"11");
        }
        else { // 4 unique symbols, unbalanced topology
            HuffNode n1 = null;
            HuffNode n2 = null;
            HuffNode n3 = null;
            HuffNode n4 = null;
            for(HuffNode l : leaves) {
                if(n1 == null || l.count > n1.count) {
                    n4 = n3;
                    n3 = n2;
                    n2 = n1;
                    n1 = l;
                }
                else if(n2 == null || l.count > n2.count) {
                    n4 = n3;
                    n3 = n2;
                    n2 = l;
                }
                else if(n3 == null || l.count > n3.count) {
                    n4 = n3;
                    n3 = l;
                }
                else {
                    n4 = l;
                }
            }
            nucToBin.put(n1.symbol,"1");
            nucToBin.put(n2.symbol,"01");
            nucToBin.put(n3.symbol,"001");
            nucToBin.put(n4.symbol,"000");
            infoByte = "00";    // 4 unique symbols
            switch(n1.symbol) { // which symbol is largest?
                case 'A': infoByte += "00"; break;
                case 'C': infoByte += "01"; break;
                case 'G': infoByte += "10"; break;
                case 'T': infoByte += "11"; break;
                default: System.exit(-1);
            }
            switch(n2.symbol) { // which symbol is 2nd largest?
                case 'A': infoByte += "00"; break;
                case 'C': infoByte += "01"; break;
                case 'G': infoByte += "10"; break;
                case 'T': infoByte += "11"; break;
                default: System.exit(-1);
            }
            switch(n3.symbol) { // which symbol is 3rd largest? smallest can be inferred from these 3
                case 'A': infoByte += "00"; break;
                case 'C': infoByte += "01"; break;
                case 'G': infoByte += "10"; break;
                case 'T': infoByte += "11"; break;
                default: System.exit(-1);
            }
        }
        
        // encode file
        try {
            InputStream in = new BufferedInputStream(new FileInputStream(new File(INFILE)));
            DataOutputStream out = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(new File(OUTFILE))));
            
            // write infoByte, then numChars
            out.writeByte(Integer.parseInt(infoByte, 2));
            out.writeLong(numChars);
            
            // if one unique symbol, we're done with just infobyte and numchars
            if(numUnique == 1) {
                in.close();
                out.close();
                System.exit(0);
            }
            
            String buf = "";
            int b;
            while((b = in.read()) != -1) {
                buf += nucToBin.get((char)b);
                if(buf.length() >= 8) {
                    out.writeByte(Integer.parseInt(buf.substring(0,8),2));
                    buf = buf.substring(8);
                }
            }
            if(buf.length() > 0) {
                while(buf.length() < 8) {
                    buf += "0";
                }
                out.writeByte(Integer.parseInt(buf,2));
            }
            out.close();
        }
        catch(FileNotFoundException e) {
            System.err.println("ERROR: File \"" + INFILE + "\" not found!");
            System.exit(-1);
        }
        catch(IOException e) {
            System.err.println("ERROR: IOException while reading \"" + INFILE + "\"!");
            System.exit(-1);
        }
    }
    
    /* HUFFMAN DECOMPRESS */
    public static void decompress(String INFILE, String OUTFILE) {
        // read file and decode
        try {
            DataInputStream in = new DataInputStream(new BufferedInputStream(new FileInputStream(new File(INFILE))));
            DataOutputStream out = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(new File(OUTFILE))));
            byte infoByte = in.readByte();
            int numUnique = (int)(infoByte >>> 6) & 3;
            
            // if 1 unique symbol, just read length and print that symbol length times
            if(numUnique == 1) {
                char symbol = 'Z';
                switch((int)(infoByte >>> 4) & 3) { // the symbol
                    case 0: symbol = 'A'; break;
                    case 1: symbol = 'C'; break;
                    case 2: symbol = 'G'; break;
                    case 3: symbol = 'T'; break;
                    default: System.exit(-1);
                }
                long numChars = in.readLong();
                for(long i = 0; i < numChars; ++i) {
                    out.writeByte((byte)symbol);
                }
                in.close();
                out.close();
                System.exit(0);
            }
            
            // build Huffman Tree
            HuffNode root = new HuffNode('\0');
            if((int)infoByte == 0) { // all 4 symbols, topology 1, so balanced
                root.l = new HuffNode('\0');
                root.r = new HuffNode('\0');
                root.r.r = new HuffNode('A'); // 00 = A
                root.r.l = new HuffNode('C'); // 01 = C
                root.l.r = new HuffNode('G'); // 10 = G
                root.l.l = new HuffNode('T'); // 11 = T
            }
            else if(numUnique == 2) { // 2 unique symbols
                switch((int)(infoByte >>> 4) & 3) { // the larger symbol
                    case 0: root.l = new HuffNode('A'); break;
                    case 1: root.l = new HuffNode('C'); break;
                    case 2: root.l = new HuffNode('G'); break;
                    case 3: root.l = new HuffNode('T'); break;
                    default: System.exit(-1);
                }
                switch((int)(infoByte >>> 2) & 3) { // the smaller symbol
                    case 0: root.r = new HuffNode('A'); break;
                    case 1: root.r = new HuffNode('C'); break;
                    case 2: root.r = new HuffNode('G'); break;
                    case 3: root.r = new HuffNode('T'); break;
                    default: System.exit(-1);
                }
            }
            else if(numUnique == 3) { // 3 unique symbols
                root.r = new HuffNode('\0');
                switch((int)(infoByte >>> 4) & 3) { // the largest symbol
                    case 0: root.l = new HuffNode('A'); break;
                    case 1: root.l = new HuffNode('C'); break;
                    case 2: root.l = new HuffNode('G'); break;
                    case 3: root.l = new HuffNode('T'); break;
                    default: System.exit(-1);
                }
                switch((int)(infoByte >>> 2) & 3) { // the middle symbol
                    case 0: root.r.l = new HuffNode('A'); break;
                    case 1: root.r.l = new HuffNode('C'); break;
                    case 2: root.r.l = new HuffNode('G'); break;
                    case 3: root.r.l = new HuffNode('T'); break;
                    default: System.exit(-1);
                }
                switch((int)infoByte & 3) { // the smallest symbol
                    case 0: root.r.r = new HuffNode('A'); break;
                    case 1: root.r.r = new HuffNode('C'); break;
                    case 2: root.r.r = new HuffNode('G'); break;
                    case 3: root.r.r = new HuffNode('T'); break;
                    default: System.exit(-1);
                }
            }
            else { // 4 unique symbols, unbalanced topology
                char l1 = 'Z';
                char l2 = 'Z';
                char l3 = 'Z';
                boolean[] used = new boolean[4];
                switch((int)(infoByte >>> 4) & 3) { // extract l1
                    case 0: l1 = 'A'; used[0] = true; break;
                    case 1: l1 = 'C'; used[1] = true; break;
                    case 2: l1 = 'G'; used[2] = true; break;
                    case 3: l1 = 'T'; used[3] = true; break;
                    default: System.exit(-1);
                }
                switch((int)(infoByte >>> 2) & 3) { // extract l2
                    case 0: l2 = 'A'; used[0] = true; break;
                    case 1: l2 = 'C'; used[1] = true; break;
                    case 2: l2 = 'G'; used[2] = true; break;
                    case 3: l2 = 'T'; used[3] = true; break;
                    default: System.exit(-1);
                }
                switch((int)infoByte & 3) { // extract l3
                    case 0: l3 = 'A'; used[0] = true; break;
                    case 1: l3 = 'C'; used[1] = true; break;
                    case 2: l3 = 'G'; used[2] = true; break;
                    case 3: l3 = 'T'; used[3] = true; break;
                    default: System.exit(-1);
                }
                root.l = new HuffNode(l1); // length-1 symbol
                root.r = new HuffNode('\0');
                root.r.l = new HuffNode(l2); // length-2 symbol
                root.r.r = new HuffNode('\0');
                root.r.r.l = new HuffNode(l3); // length-3 symbol
                for(int i = 0; i < 4; ++i) {
                    if(!used[i]) {
                        switch(i) {
                            case 0: root.r.r.r = new HuffNode('A'); break;
                            case 1: root.r.r.r = new HuffNode('C'); break;
                            case 2: root.r.r.r = new HuffNode('G'); break;
                            case 3: root.r.r.r = new HuffNode('T'); break;
                        }
                    }
                }
            }
            
            // read rest of file and decode
            long numChars = in.readLong();
            long printed = 0L;
            HuffNode c = root;
            while(true) {
                byte buf = in.readByte();
                for(int i = 7; i >= 0; --i) {
                    int bit = (buf >>> i) & 1;
                    if(bit == 0) {
                        if(c.r != null) {
                            c = c.r;
                            if(c.l == null && c.r == null) {
                                out.writeByte((byte)(c.symbol));
                                if(++printed == numChars) {
                                    in.close();
                                    out.close();
                                    System.exit(0);
                                }
                                c = root;
                            }
                        }
                        else {
                            c = root.r;
                        }
                    }
                    else {
                        if(c.l != null) {
                            c = c.l;
                            if(c.l == null && c.r == null) {
                                out.writeByte((byte)(c.symbol));
                                if(++printed == numChars) {
                                    in.close();
                                    out.close();
                                    System.exit(0);
                                }
                                c = root;
                            }
                        }
                        else {
                            c = root.l;
                        }
                    }
                }
            }
        }
        catch(Exception e) {
            e.printStackTrace();
            System.exit(-1);
        }
    }
}

/* HELPER CLASS: HuffNode */
class HuffNode implements Comparable<HuffNode>{
    public char symbol;
    public long count;
    public HuffNode l; // 1 child
    public HuffNode r; // 0 child
    public HuffNode p; // parent
    public String code;
    
    public HuffNode(char s) {
        symbol = s;
        count = 0L;
        l = null;
        r = null;
        p = null;
        code = "";
    }
    
    public String code() {
        if(code != "") {
            return code;
        }
        HuffNode curr = this;
        while(curr.p != null) {
            if(curr == curr.p.l) {
                code = "1" + code;
            }
            else if(curr == curr.p.r) {
                code = "0" + code;
            }
            else {
                System.err.println("ERROR: Node is not left nor right child of parent!");
                System.exit(-1);
            }
            curr = curr.p;
        }
        return code;
    }
    
    public void print() {
        if(l != null) {
            System.out.println(this + " -1-> " + l);
            l.print();
        }
        if(r != null) {
            System.out.println(this + " -0-> " + r);
            r.print();
        }
    }
    
    @Override
    public int compareTo(HuffNode other) {
        return new Long(this.count).compareTo(new Long(other.count));
    }
    
    @Override
    public String toString() {
        return "(" + symbol + "," + count + ")";
    }
}