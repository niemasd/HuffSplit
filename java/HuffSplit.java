/* AUTHOR: Niema Moshiri
 * DNA Split Huffman Compression
 *
 * USAGE:
 * -Compress:   java HuffSplit compress <in_file>
 * -Decompress: java HuffSplit decompress <huffsplit_file>
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
import java.io.*;
import java.nio.file.Paths;
import java.nio.file.Files;
import java.util.*;

public class HuffSplit {
    // instance variables
    public static final int NUMTOPS = 165;          // number of possible topologies
    @SuppressWarnings("unchecked")
    public static HashMap<Character,String>[] TOPS = new HashMap[NUMTOPS]; // topologies
    
    /* Main Method
     */
    public static void main( String[] args ) {
        // parse arguments
        if(args.length != 2) {
            System.err.println("ERROR: Incorrect number of arguments");
            System.err.println("See file header for usage information");
            System.exit(-1);
        }
        final String IN = args[1];
        
        // get tree topologies
        for(int i = 0; i < NUMTOPS; ++i) {
            TOPS[i] = getCode(i);
        }
        
        // run relevant function
        switch(args[0]) {
            case "compress": compress(IN,IN+".hsf"); break;
            case "decompress": decompress(IN,IN.substring(0,IN.lastIndexOf('.'))); break;
            default: System.err.println("ERROR: First argument must be \"compress\" or \"decompress\"!"); System.err.println("See file header for usage information"); System.exit(-1);
        }
    }
    
    /* Compress the input file using my split Huffman algorithm
     * INPUT:  A DNA string to compress
     * OUTPUT: The compressed results of my split Huffman algorithm
     */
    public static void compress( String INFILE, String OUTFILE ) {
        // read input file
        String in = null;
        try {
            in = new String(Files.readAllBytes(Paths.get(INFILE)));
        } catch(FileNotFoundException e) {
            System.err.println("ERROR: File \"" + INFILE + "\" not found!"); System.exit(-1);
        } catch(IOException e) {
            System.err.println("ERROR: IOException while reading \"" + INFILE + "\"!"); System.exit(-1);
        }
        if(in == null) {
            System.err.println("ERROR: Something went wrong while reading \"" + INFILE + "\"!"); System.exit(-1);
        }
        if(in.length() == 0) {
            System.err.println("ERROR: Empty file!"); System.exit(-1);
        }
        if(in.charAt(0) != 'A' && in.charAt(0) != 'C' && in.charAt(0) != 'G' && in.charAt(0) != 'T' && in.charAt(0) != 'N') {
            System.err.println("ERROR: Invalid symbol: " + in.charAt(0)); System.exit(-1);
        }
        final int L = in.length();
        
        // get optimal cuts
        int[][] C = new int[2][NUMTOPS];
        byte[][] backtrack = new byte[L][NUMTOPS];
        int[] bestT = {-1,-1};
        int[] bestC = {-1,-1};
        for(int t = 0; t < NUMTOPS; ++t) {
            if(TOPS[t].containsKey(in.charAt(0))) {
                C[0][t] = 72 + TOPS[t].get(in.charAt(0)).length();
                backtrack[0][t] = (byte)t;
                if(bestC[0] == -1 || bestC[0] > C[0][t]) {
                    bestC[0] = C[0][t];
                    bestT[0] = t;
                }
            }
            else {
                C[0][t] = -1;
                backtrack[0][t] = (byte)-1;
            }
        }
        for(int i = 1; i < L; ++i) {
            int i1 = i%2;
            int i0 = (i-1)%2;
            char c = in.charAt(i);
            if(c != 'A' && c != 'C' && c != 'G' && c != 'T' && c != 'N') {
                System.err.println("ERROR: Invalid symbol: " + c); System.exit(-1);
            }
            for(int top = 0; top < NUMTOPS; ++top) {
                if(TOPS[top].containsKey(c)) {
                    int bits = TOPS[top].get(c).length();
                    if(bestT[i0] == top) {
                        C[i1][top] = C[i0][top] + bits;
                        backtrack[i][top] = (byte)top;
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
                            backtrack[i][top] = (byte)bestT[i0];
                        }
                        else {
                            C[i1][top] = sameC;
                            backtrack[i][top] = (byte)top;
                        }
                    }
                    if(bestC[i1] == -1 || bestC[i1] > C[i1][top]) {
                        bestC[i1] = C[i1][top];
                        bestT[i1] = top;
                    }
                }
                else {
                    C[i1][top] = -1;
                    backtrack[i][top] = (byte)-1;
                }
            }
            bestC[i0] = -1;
            bestT[i0] = -1;
        }
        
        // reconstruct topology path from backtrack
        int[] path = new int[L]; // path[i] is the topology we are in at character i of the input
        int LAST = (L-1)%2;
        for(int t = 0; t < NUMTOPS; ++t) {
            if(C[LAST][t] != -1 && (C[LAST][path[L-1]] == -1 || C[LAST][t] < C[LAST][path[L-1]])) {
                path[L-1] = t;
            }
        }
        for(int i = L-2; i >= 0; --i) {
            path[i] = (int)(backtrack[i+1][path[i+1]] & 0xFF);
        }
        C = null; // don't need anymore, so save memory
        ArrayList<Integer> cuts = new ArrayList<Integer>();
        cuts.add(0);
        for(int i = 1; i < L; ++i) {
            if(path[i] != path[i-1]) {
                cuts.add(i);
            }
        }
        cuts.add(L);
        
        // encode file
        try {
            DataOutputStream out = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(new File(OUTFILE))));
            for(int cut = 1; cut < cuts.size(); ++cut) {
                int start = cuts.get(cut-1);
                int end = cuts.get(cut);
                out.writeByte(path[start]); // infobyte (topology)
                out.writeInt(end-start);    // numChars
                
                // if only 1 unique symbol, only need first 5 bytes
                if(path[start] < 5) {
                    continue;
                }
                
                // encode substring
                String buf = "";
                for(int i = start; i < end; ++i) {
                    buf += TOPS[path[start]].get(in.charAt(i));
                    if(buf.length() >= 8) {
                        out.writeByte(Integer.parseInt(buf.substring(0,8),2));
                        buf = buf.substring(8);
                    }
                }
                
                // clear buffer
                if(buf.length() > 0) {
                    while(buf.length() < 8) {
                        buf += "0";
                    }
                    out.writeByte(Integer.parseInt(buf,2));
                }
            }
            out.close();
        } catch(Exception e) {
            e.printStackTrace();
            System.exit(-1);
        }
    }
    
    /* Decompress the input files (regular Huffman decompression on each)
     * INPUT:  The prefix of the files to decompress
     * OUTPUT: The uncompressed file
     */
    public static void decompress( String INFILE, String OUTFILE ) {
        DataInputStream in = null;
        DataOutputStream out = null;
        try {
            // set up files
            in = new DataInputStream(new BufferedInputStream(new FileInputStream(new File(INFILE))));
            out = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(new File(OUTFILE))));
            
            // decompress file
            while(true) {
                int top = (in.readByte() & 0xFF);
                int numChars = in.readInt();
                if(top < 5) {
                    char symbol = 'Z';
                    switch(top) {
                        case 0: symbol = 'A'; break;
                        case 1: symbol = 'C'; break;
                        case 2: symbol = 'G'; break;
                        case 3: symbol = 'T'; break;
                        case 4: symbol = 'N'; break;
                        default: System.err.println("ERROR: Unrecognized topology: " + top); System.exit(-1);
                    }
                    for(int i = 0; i < numChars; ++i) {
                        out.writeByte((byte)symbol);
                        continue;
                    }
                }
                else {
                    Node root = buildTree(top);
                    Node c = root;
                    int printed = 0;
                    while(printed < numChars) {
                        byte buf = in.readByte();
                        for(int i = 7; i >= 0; --i) {
                            int bit = (buf >>> i) & 1;
                            if(bit == 0) {
                                if(c.c0 != null) {
                                    c = c.c0;
                                    if(c.c1 == null && c.c0 == null) {
                                        out.writeByte((byte)(c.symbol));
                                        if(++printed == numChars) {
                                            break;
                                        }
                                        c = root;
                                    }
                                }
                                else {
                                    System.err.println("ERROR: Invalid c0");
                                    System.exit(-1);
                                }
                            }
                            else {
                                if(c.c1 != null) {
                                    c = c.c1;
                                    if(c.c1 == null && c.c0 == null) {
                                        out.writeByte((byte)(c.symbol));
                                        if(++printed == numChars) {
                                            break;
                                        }
                                        c = root;
                                    }
                                }
                                else {
                                    System.err.println("ERROR: Invalid c1");
                                    System.exit(-1);
                                }
                            }
                        }
                    }
                }
            }
        } catch(EOFException e) {
        } catch(FileNotFoundException e) {
            System.err.println("ERROR: One of the files was not found!"); e.printStackTrace(); System.exit(-1);
        } catch(IOException e) {
            System.err.println("ERROR: IO Exception!"); e.printStackTrace(); System.exit(-1);
        }
        try {
            if(in != null) {
                in.close();
                out.close();
            }
        } catch(IOException e) {
            System.err.println("ERROR: IO Exception when closing input or output file!"); e.printStackTrace(); System.exit(-1);
        }
    }
    
    /* Given an integer, return the corresponding tree (see Topologies.pptx)
     * INPUT:  An integer (between 0 and 164, inclusive)
     * OUTPUT: The root node of the corresponding tree
     */
    public static Node buildTree( int topology ) {
        HashMap<Character,String> code = getCode(topology);
        Node root = new Node((char)0);
        if(topology < 5) {
            return root;
        }
        else if(topology < 15) {
            for(Character c : code.keySet()) {
                if(code.get(c).equals("1")) {
                    root.c1 = new Node(c);
                }
                else {
                    root.c0 = new Node(c);
                }
            }
        }
        else if(topology < 45) {
            Node nodeA = null;
            Node nodeB = null;
            Node nodeC = null;
            for(Character c : code.keySet()) {
                if(code.get(c).equals("1")) {
                    nodeA = new Node(c);
                }
                else if(code.get(c).equals("01")) {
                    nodeB = new Node(c);
                }
                else {
                    nodeC = new Node(c);
                }
            }
            root.c1 = nodeA;
            root.c0 = new Node((char)0);
            root.c0.c1 = nodeB;
            root.c0.c0 = nodeC;
        }
        else if(topology < 50) {
            Node nodeA = null;
            Node nodeB = null;
            Node nodeC = null;
            Node nodeD = null;
            for(Character c : code.keySet()) {
                if(code.get(c).equals("11")) {
                    nodeA = new Node(c);
                }
                else if(code.get(c).equals("10")) {
                    nodeB = new Node(c);
                }
                else if(code.get(c).equals("01")) {
                    nodeC = new Node(c);
                }
                else {
                    nodeD = new Node(c);
                }
            }
            root.c1 = new Node((char)0);
            root.c0 = new Node((char)0);
            root.c1.c1 = nodeA;
            root.c1.c0 = nodeB;
            root.c0.c1 = nodeC;
            root.c0.c0 = nodeD;
        }
        else if(topology < 90) {
            Node nodeA = null;
            Node nodeB = null;
            Node nodeC = null;
            Node nodeD = null;
            for(Character c : code.keySet()) {
                if(code.get(c).equals("1")) {
                    nodeA = new Node(c);
                }
                else if(code.get(c).equals("01")) {
                    nodeB = new Node(c);
                }
                else if(code.get(c).equals("001")) {
                    nodeC = new Node(c);
                }
                else {
                    nodeD = new Node(c);
                }
            }
            root.c1 = nodeA;
            root.c0 = new Node((char)0);
            root.c0.c1 = nodeB;
            root.c0.c0 = new Node((char)0);
            root.c0.c0.c1 = nodeC;
            root.c0.c0.c0 = nodeD;
        }
        else if(topology < 150) {
            Node nodeA = null;
            Node nodeB = null;
            Node nodeC = null;
            Node nodeD = null;
            Node nodeE = null;
            for(Character c : code.keySet()) {
                if(code.get(c).equals("1")) {
                    nodeA = new Node(c);
                }
                else if(code.get(c).equals("01")) {
                    nodeB = new Node(c);
                }
                else if(code.get(c).equals("001")) {
                    nodeC = new Node(c);
                }
                else if(code.get(c).equals("0001")) {
                    nodeD = new Node(c);
                }
                else {
                    nodeE = new Node(c);
                }
            }
            root.c1 = nodeA;
            root.c0 = new Node((char)0);
            root.c0.c1 = nodeB;
            root.c0.c0 = new Node((char)0);
            root.c0.c0.c1 = nodeC;
            root.c0.c0.c0 = new Node((char)0);
            root.c0.c0.c0.c1 = nodeD;
            root.c0.c0.c0.c0 = nodeE;
        }
        else if(topology < 160) {
            Node nodeA = null;
            Node nodeB = null;
            Node nodeC = null;
            Node nodeD = null;
            Node nodeE = null;
            for(Character c : code.keySet()) {
                if(code.get(c).equals("11")) {
                    nodeA = new Node(c);
                }
                else if(code.get(c).equals("10")) {
                    nodeB = new Node(c);
                }
                else if(code.get(c).equals("01")) {
                    nodeC = new Node(c);
                }
                else if(code.get(c).equals("001")) {
                    nodeD = new Node(c);
                }
                else {
                    nodeE = new Node(c);
                }
            }
            root.c1 = new Node((char)0);
            root.c1.c1 = nodeA;
            root.c1.c0 = nodeB;
            root.c0 = new Node((char)0);
            root.c0.c1 = nodeC;
            root.c0.c0 = new Node((char)0);
            root.c0.c0.c1 = nodeD;
            root.c0.c0.c0 = nodeE;
        }
        else if(topology < 165) {
            Node nodeA = null;
            Node nodeB = null;
            Node nodeC = null;
            Node nodeD = null;
            Node nodeE = null;
            for(Character c : code.keySet()) {
                if(code.get(c).equals("1")) {
                    nodeA = new Node(c);
                }
                else if(code.get(c).equals("011")) {
                    nodeB = new Node(c);
                }
                else if(code.get(c).equals("010")) {
                    nodeC = new Node(c);
                }
                else if(code.get(c).equals("001")) {
                    nodeD = new Node(c);
                }
                else {
                    nodeE = new Node(c);
                }
            }
            root.c1 = nodeA;
            root.c0 = new Node((char)0);
            root.c0.c1 = new Node((char)0);
            root.c0.c1.c1 = nodeB;
            root.c0.c1.c0 = nodeC;
            root.c0.c0 = new Node((char)0);
            root.c0.c0.c1 = nodeD;
            root.c0.c0.c0 = nodeE;
        }
        else {
            System.err.println("ERROR: Trying to get tree for topology > 164!");
            System.exit(-1);
        }
        return root;
    }
    
    /* Given an integer, return the corresponding topology (see Topologies.pptx)
     * INPUT:  An integer (between 0 and 164, inclusive)
     * OUTPUT: A HashMap containing the code of the corresponding topology
     */
    public static HashMap<Character,String> getCode( int topology ) {
        HashMap<Character,String> code = new HashMap<Character,String>();
        switch(topology) {
            // 1 Unique Character
            case 0:   code.put('A',""); break;
            case 1:   code.put('C',""); break;
            case 2:   code.put('G',""); break;
            case 3:   code.put('T',""); break;
            case 4:   code.put('N',""); break;
            
            // 2 Unique Characters
            case 5:   code.put('A',"1"); code.put('C',"0"); break;
            case 6:   code.put('A',"1"); code.put('G',"0"); break;
            case 7:   code.put('A',"1"); code.put('T',"0"); break;
            case 8:   code.put('A',"1"); code.put('N',"0"); break;
            case 9:   code.put('C',"1"); code.put('G',"0"); break;
            case 10:  code.put('C',"1"); code.put('T',"0"); break;
            case 11:  code.put('C',"1"); code.put('N',"0"); break;
            case 12:  code.put('G',"1"); code.put('T',"0"); break;
            case 13:  code.put('G',"1"); code.put('N',"0"); break;
            case 14:  code.put('T',"1"); code.put('N',"0"); break;
            
            // 3 Unique Characters
            case 15:  code.put('A',"1"); code.put('C',"01"); code.put('G',"00"); break;
            case 16:  code.put('A',"1"); code.put('C',"01"); code.put('T',"00"); break;
            case 17:  code.put('A',"1"); code.put('C',"01"); code.put('N',"00"); break;
            case 18:  code.put('A',"1"); code.put('G',"01"); code.put('T',"00"); break;
            case 19:  code.put('A',"1"); code.put('G',"01"); code.put('N',"00"); break;
            case 20:  code.put('A',"1"); code.put('T',"01"); code.put('N',"00"); break;
            case 21:  code.put('C',"1"); code.put('A',"01"); code.put('G',"00"); break;
            case 22:  code.put('C',"1"); code.put('A',"01"); code.put('T',"00"); break;
            case 23:  code.put('C',"1"); code.put('A',"01"); code.put('N',"00"); break;
            case 24:  code.put('C',"1"); code.put('G',"01"); code.put('T',"00"); break;
            case 25:  code.put('C',"1"); code.put('G',"01"); code.put('N',"00"); break;
            case 26:  code.put('C',"1"); code.put('T',"01"); code.put('N',"00"); break;
            case 27:  code.put('G',"1"); code.put('A',"01"); code.put('C',"00"); break;
            case 28:  code.put('G',"1"); code.put('A',"01"); code.put('T',"00"); break;
            case 29:  code.put('G',"1"); code.put('A',"01"); code.put('N',"00"); break;
            case 30:  code.put('G',"1"); code.put('C',"01"); code.put('T',"00"); break;
            case 31:  code.put('G',"1"); code.put('C',"01"); code.put('N',"00"); break;
            case 32:  code.put('G',"1"); code.put('T',"01"); code.put('N',"00"); break;
            case 33:  code.put('T',"1"); code.put('A',"01"); code.put('C',"00"); break;
            case 34:  code.put('T',"1"); code.put('A',"01"); code.put('G',"00"); break;
            case 35:  code.put('T',"1"); code.put('A',"01"); code.put('N',"00"); break;
            case 36:  code.put('T',"1"); code.put('C',"01"); code.put('G',"00"); break;
            case 37:  code.put('T',"1"); code.put('C',"01"); code.put('N',"00"); break;
            case 38:  code.put('T',"1"); code.put('G',"01"); code.put('N',"00"); break;
            case 39:  code.put('N',"1"); code.put('A',"01"); code.put('C',"00"); break;
            case 40:  code.put('N',"1"); code.put('A',"01"); code.put('G',"00"); break;
            case 41:  code.put('N',"1"); code.put('A',"01"); code.put('T',"00"); break;
            case 42:  code.put('N',"1"); code.put('C',"01"); code.put('G',"00"); break;
            case 43:  code.put('N',"1"); code.put('C',"01"); code.put('T',"00"); break;
            case 44:  code.put('N',"1"); code.put('G',"01"); code.put('T',"00"); break;
            
            // 4 Unique Characters (Balanced)
            case 45:  code.put('A',"11"); code.put('C',"10"); code.put('G',"01"); code.put('T',"00"); break;
            case 46:  code.put('A',"11"); code.put('C',"10"); code.put('G',"01"); code.put('N',"00"); break;
            case 47:  code.put('A',"11"); code.put('C',"10"); code.put('T',"01"); code.put('N',"00"); break;
            case 48:  code.put('A',"11"); code.put('G',"10"); code.put('T',"01"); code.put('N',"00"); break;
            case 49:  code.put('C',"11"); code.put('G',"10"); code.put('T',"01"); code.put('N',"00"); break;
            
            // 4 Unique Characters (Unbalanced)
            case 50:  code.put('A',"1"); code.put('C',"01"); code.put('G',"001"); code.put('T',"000"); break;
            case 51:  code.put('A',"1"); code.put('C',"01"); code.put('G',"001"); code.put('N',"000"); break;
            case 52:  code.put('A',"1"); code.put('G',"01"); code.put('C',"001"); code.put('T',"000"); break;
            case 53:  code.put('A',"1"); code.put('G',"01"); code.put('C',"001"); code.put('N',"000"); break;
            case 54:  code.put('A',"1"); code.put('T',"01"); code.put('C',"001"); code.put('G',"000"); break;
            case 55:  code.put('A',"1"); code.put('T',"01"); code.put('C',"001"); code.put('N',"000"); break;
            case 56:  code.put('A',"1"); code.put('N',"01"); code.put('C',"001"); code.put('G',"000"); break;
            case 57:  code.put('A',"1"); code.put('N',"01"); code.put('C',"001"); code.put('T',"000"); break;
            case 58:  code.put('C',"1"); code.put('A',"01"); code.put('G',"001"); code.put('T',"000"); break;
            case 59:  code.put('C',"1"); code.put('A',"01"); code.put('G',"001"); code.put('N',"000"); break;
            case 60:  code.put('C',"1"); code.put('G',"01"); code.put('A',"001"); code.put('T',"000"); break;
            case 61:  code.put('C',"1"); code.put('G',"01"); code.put('A',"001"); code.put('N',"000"); break;
            case 62:  code.put('C',"1"); code.put('T',"01"); code.put('A',"001"); code.put('G',"000"); break;
            case 63:  code.put('C',"1"); code.put('T',"01"); code.put('A',"001"); code.put('N',"000"); break;
            case 64:  code.put('C',"1"); code.put('N',"01"); code.put('A',"001"); code.put('G',"000"); break;
            case 65:  code.put('C',"1"); code.put('N',"01"); code.put('A',"001"); code.put('T',"000"); break;
            case 66:  code.put('G',"1"); code.put('A',"01"); code.put('C',"001"); code.put('T',"000"); break;
            case 67:  code.put('G',"1"); code.put('A',"01"); code.put('C',"001"); code.put('N',"000"); break;
            case 68:  code.put('G',"1"); code.put('C',"01"); code.put('A',"001"); code.put('T',"000"); break;
            case 69:  code.put('G',"1"); code.put('C',"01"); code.put('A',"001"); code.put('N',"000"); break;
            case 70:  code.put('G',"1"); code.put('T',"01"); code.put('A',"001"); code.put('C',"000"); break;
            case 71:  code.put('G',"1"); code.put('T',"01"); code.put('A',"001"); code.put('N',"000"); break;
            case 72:  code.put('G',"1"); code.put('N',"01"); code.put('A',"001"); code.put('C',"000"); break;
            case 73:  code.put('G',"1"); code.put('N',"01"); code.put('A',"001"); code.put('T',"000"); break;
            case 74:  code.put('T',"1"); code.put('A',"01"); code.put('C',"001"); code.put('G',"000"); break;
            case 75:  code.put('T',"1"); code.put('A',"01"); code.put('C',"001"); code.put('N',"000"); break;
            case 76:  code.put('T',"1"); code.put('C',"01"); code.put('A',"001"); code.put('G',"000"); break;
            case 77:  code.put('T',"1"); code.put('C',"01"); code.put('A',"001"); code.put('N',"000"); break;
            case 78:  code.put('T',"1"); code.put('G',"01"); code.put('A',"001"); code.put('C',"000"); break;
            case 79:  code.put('T',"1"); code.put('G',"01"); code.put('A',"001"); code.put('N',"000"); break;
            case 80:  code.put('T',"1"); code.put('N',"01"); code.put('A',"001"); code.put('C',"000"); break;
            case 81:  code.put('T',"1"); code.put('N',"01"); code.put('A',"001"); code.put('G',"000"); break;
            case 82:  code.put('N',"1"); code.put('A',"01"); code.put('C',"001"); code.put('G',"000"); break;
            case 83:  code.put('N',"1"); code.put('A',"01"); code.put('C',"001"); code.put('T',"000"); break;
            case 84:  code.put('N',"1"); code.put('C',"01"); code.put('A',"001"); code.put('G',"000"); break;
            case 85:  code.put('N',"1"); code.put('C',"01"); code.put('A',"001"); code.put('T',"000"); break;
            case 86:  code.put('N',"1"); code.put('G',"01"); code.put('A',"001"); code.put('C',"000"); break;
            case 87:  code.put('N',"1"); code.put('G',"01"); code.put('A',"001"); code.put('T',"000"); break;
            case 88:  code.put('N',"1"); code.put('T',"01"); code.put('A',"001"); code.put('C',"000"); break;
            case 89:  code.put('N',"1"); code.put('T',"01"); code.put('A',"001"); code.put('G',"000"); break;
            
            // 5 Unique Characters (Line)
            case 90:  code.put('A',"1"); code.put('C',"01"); code.put('G',"001"); code.put('T',"0001"); code.put('N',"0000"); break;
            case 91:  code.put('A',"1"); code.put('C',"01"); code.put('T',"001"); code.put('G',"0001"); code.put('N',"0000"); break;
            case 92:  code.put('A',"1"); code.put('C',"01"); code.put('N',"001"); code.put('G',"0001"); code.put('T',"0000"); break;
            case 93:  code.put('A',"1"); code.put('G',"01"); code.put('C',"001"); code.put('T',"0001"); code.put('N',"0000"); break;
            case 94:  code.put('A',"1"); code.put('G',"01"); code.put('T',"001"); code.put('C',"0001"); code.put('N',"0000"); break;
            case 95:  code.put('A',"1"); code.put('G',"01"); code.put('N',"001"); code.put('C',"0001"); code.put('T',"0000"); break;
            case 96:  code.put('A',"1"); code.put('T',"01"); code.put('C',"001"); code.put('G',"0001"); code.put('N',"0000"); break;
            case 97:  code.put('A',"1"); code.put('T',"01"); code.put('G',"001"); code.put('C',"0001"); code.put('N',"0000"); break;
            case 98:  code.put('A',"1"); code.put('T',"01"); code.put('N',"001"); code.put('C',"0001"); code.put('G',"0000"); break;
            case 99:  code.put('A',"1"); code.put('N',"01"); code.put('C',"001"); code.put('G',"0001"); code.put('T',"0000"); break;
            case 100: code.put('A',"1"); code.put('N',"01"); code.put('G',"001"); code.put('C',"0001"); code.put('T',"0000"); break;
            case 101: code.put('A',"1"); code.put('N',"01"); code.put('T',"001"); code.put('C',"0001"); code.put('G',"0000"); break;
            case 102: code.put('C',"1"); code.put('A',"01"); code.put('G',"001"); code.put('T',"0001"); code.put('N',"0000"); break;
            case 103: code.put('C',"1"); code.put('A',"01"); code.put('T',"001"); code.put('G',"0001"); code.put('N',"0000"); break;
            case 104: code.put('C',"1"); code.put('A',"01"); code.put('N',"001"); code.put('G',"0001"); code.put('T',"0000"); break;
            case 105: code.put('C',"1"); code.put('G',"01"); code.put('A',"001"); code.put('T',"0001"); code.put('N',"0000"); break;
            case 106: code.put('C',"1"); code.put('G',"01"); code.put('T',"001"); code.put('A',"0001"); code.put('N',"0000"); break;
            case 107: code.put('C',"1"); code.put('G',"01"); code.put('N',"001"); code.put('A',"0001"); code.put('T',"0000"); break;
            case 108: code.put('C',"1"); code.put('T',"01"); code.put('A',"001"); code.put('G',"0001"); code.put('N',"0000"); break;
            case 109: code.put('C',"1"); code.put('T',"01"); code.put('G',"001"); code.put('A',"0001"); code.put('N',"0000"); break;
            case 110: code.put('C',"1"); code.put('T',"01"); code.put('N',"001"); code.put('A',"0001"); code.put('G',"0000"); break;
            case 111: code.put('C',"1"); code.put('N',"01"); code.put('A',"001"); code.put('G',"0001"); code.put('T',"0000"); break;
            case 112: code.put('C',"1"); code.put('N',"01"); code.put('G',"001"); code.put('A',"0001"); code.put('T',"0000"); break;
            case 113: code.put('C',"1"); code.put('N',"01"); code.put('T',"001"); code.put('A',"0001"); code.put('G',"0000"); break;
            case 114: code.put('G',"1"); code.put('A',"01"); code.put('C',"001"); code.put('T',"0001"); code.put('N',"0000"); break;
            case 115: code.put('G',"1"); code.put('A',"01"); code.put('T',"001"); code.put('C',"0001"); code.put('N',"0000"); break;
            case 116: code.put('G',"1"); code.put('A',"01"); code.put('N',"001"); code.put('C',"0001"); code.put('T',"0000"); break;
            case 117: code.put('G',"1"); code.put('C',"01"); code.put('A',"001"); code.put('T',"0001"); code.put('N',"0000"); break;
            case 118: code.put('G',"1"); code.put('C',"01"); code.put('T',"001"); code.put('A',"0001"); code.put('N',"0000"); break;
            case 119: code.put('G',"1"); code.put('C',"01"); code.put('N',"001"); code.put('A',"0001"); code.put('T',"0000"); break;
            case 120: code.put('G',"1"); code.put('T',"01"); code.put('A',"001"); code.put('C',"0001"); code.put('N',"0000"); break;
            case 121: code.put('G',"1"); code.put('T',"01"); code.put('C',"001"); code.put('A',"0001"); code.put('N',"0000"); break;
            case 122: code.put('G',"1"); code.put('T',"01"); code.put('N',"001"); code.put('A',"0001"); code.put('C',"0000"); break;
            case 123: code.put('G',"1"); code.put('N',"01"); code.put('A',"001"); code.put('C',"0001"); code.put('T',"0000"); break;
            case 124: code.put('G',"1"); code.put('N',"01"); code.put('C',"001"); code.put('A',"0001"); code.put('T',"0000"); break;
            case 125: code.put('G',"1"); code.put('N',"01"); code.put('T',"001"); code.put('A',"0001"); code.put('C',"0000"); break;
            case 126: code.put('T',"1"); code.put('A',"01"); code.put('C',"001"); code.put('G',"0001"); code.put('N',"0000"); break;
            case 127: code.put('T',"1"); code.put('A',"01"); code.put('G',"001"); code.put('C',"0001"); code.put('N',"0000"); break;
            case 128: code.put('T',"1"); code.put('A',"01"); code.put('N',"001"); code.put('C',"0001"); code.put('G',"0000"); break;
            case 129: code.put('T',"1"); code.put('C',"01"); code.put('A',"001"); code.put('G',"0001"); code.put('N',"0000"); break;
            case 130: code.put('T',"1"); code.put('C',"01"); code.put('G',"001"); code.put('A',"0001"); code.put('N',"0000"); break;
            case 131: code.put('T',"1"); code.put('C',"01"); code.put('N',"001"); code.put('A',"0001"); code.put('G',"0000"); break;
            case 132: code.put('T',"1"); code.put('G',"01"); code.put('A',"001"); code.put('C',"0001"); code.put('N',"0000"); break;
            case 133: code.put('T',"1"); code.put('G',"01"); code.put('C',"001"); code.put('A',"0001"); code.put('N',"0000"); break;
            case 134: code.put('T',"1"); code.put('G',"01"); code.put('N',"001"); code.put('A',"0001"); code.put('C',"0000"); break;
            case 135: code.put('T',"1"); code.put('N',"01"); code.put('A',"001"); code.put('C',"0001"); code.put('G',"0000"); break;
            case 136: code.put('T',"1"); code.put('N',"01"); code.put('C',"001"); code.put('A',"0001"); code.put('G',"0000"); break;
            case 137: code.put('T',"1"); code.put('N',"01"); code.put('G',"001"); code.put('A',"0001"); code.put('C',"0000"); break;
            case 138: code.put('N',"1"); code.put('A',"01"); code.put('C',"001"); code.put('G',"0001"); code.put('T',"0000"); break;
            case 139: code.put('N',"1"); code.put('A',"01"); code.put('G',"001"); code.put('C',"0001"); code.put('T',"0000"); break;
            case 140: code.put('N',"1"); code.put('A',"01"); code.put('T',"001"); code.put('C',"0001"); code.put('G',"0000"); break;
            case 141: code.put('N',"1"); code.put('C',"01"); code.put('A',"001"); code.put('G',"0001"); code.put('T',"0000"); break;
            case 142: code.put('N',"1"); code.put('C',"01"); code.put('G',"001"); code.put('A',"0001"); code.put('T',"0000"); break;
            case 143: code.put('N',"1"); code.put('C',"01"); code.put('T',"001"); code.put('A',"0001"); code.put('G',"0000"); break;
            case 144: code.put('N',"1"); code.put('G',"01"); code.put('A',"001"); code.put('C',"0001"); code.put('T',"0000"); break;
            case 145: code.put('N',"1"); code.put('G',"01"); code.put('C',"001"); code.put('A',"0001"); code.put('T',"0000"); break;
            case 146: code.put('N',"1"); code.put('G',"01"); code.put('T',"001"); code.put('A',"0001"); code.put('C',"0000"); break;
            case 147: code.put('N',"1"); code.put('T',"01"); code.put('A',"001"); code.put('C',"0001"); code.put('G',"0000"); break;
            case 148: code.put('N',"1"); code.put('T',"01"); code.put('C',"001"); code.put('A',"0001"); code.put('G',"0000"); break;
            case 149: code.put('N',"1"); code.put('T',"01"); code.put('G',"001"); code.put('A',"0001"); code.put('C',"0000"); break;
            
            // 5 Unique Characters (Bend 1)
            case 150: code.put('A',"11"); code.put('C',"10"); code.put('G',"01"); code.put('T',"001"); code.put('N',"000"); break;
            case 151: code.put('A',"11"); code.put('C',"10"); code.put('T',"01"); code.put('G',"001"); code.put('N',"000"); break;
            case 152: code.put('A',"11"); code.put('G',"10"); code.put('T',"01"); code.put('C',"001"); code.put('N',"000"); break;
            case 153: code.put('C',"11"); code.put('G',"10"); code.put('T',"01"); code.put('A',"001"); code.put('N',"000"); break;
            case 154: code.put('A',"11"); code.put('C',"10"); code.put('N',"01"); code.put('G',"001"); code.put('T',"000"); break;
            case 155: code.put('A',"11"); code.put('G',"10"); code.put('N',"01"); code.put('C',"001"); code.put('T',"000"); break;
            case 156: code.put('C',"11"); code.put('G',"10"); code.put('N',"01"); code.put('A',"001"); code.put('T',"000"); break;
            case 157: code.put('A',"11"); code.put('T',"10"); code.put('N',"01"); code.put('C',"001"); code.put('G',"000"); break;
            case 158: code.put('C',"11"); code.put('T',"10"); code.put('N',"01"); code.put('A',"001"); code.put('G',"000"); break;
            case 159: code.put('G',"11"); code.put('T',"10"); code.put('N',"01"); code.put('A',"001"); code.put('C',"000"); break;
            
            // 5 Unique Characters (Bend 2)
            case 160: code.put('A',"1"); code.put('C',"011"); code.put('G',"010"); code.put('T',"001"); code.put('N',"000"); break;
            case 161: code.put('C',"1"); code.put('A',"011"); code.put('G',"010"); code.put('T',"001"); code.put('N',"000"); break;
            case 162: code.put('G',"1"); code.put('A',"011"); code.put('C',"010"); code.put('T',"001"); code.put('N',"000"); break;
            case 163: code.put('T',"1"); code.put('A',"011"); code.put('C',"010"); code.put('G',"001"); code.put('N',"000"); break;
            case 164: code.put('N',"1"); code.put('A',"011"); code.put('C',"010"); code.put('G',"001"); code.put('T',"000"); break;
            
            // Out of Bounds
            default: System.err.println("INVALID TREE TOPOLOGY"); System.exit(-1);
        }
        return code;
    }
}

/* Helper Class: Node
 */
class Node {
    public char symbol;
    public Node c0;
    public Node c1;
    public Node(char s) {
        symbol = s; c0 = null; c1 = null;
    }
}
