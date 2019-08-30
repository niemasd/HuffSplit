# USAGE: python randomDNA.py <length> [P(A) P(C) P(G) P(T)]
import sys,random

def weightedChoice(choices):
    total = sum(w for c, w in choices)
    r = random.uniform(0,total)
    upto = 0
    for c,w in choices:
        if upto + w >= r:
            return c
        upto += w

if len(sys.argv) == 2:
    sys.stdout.write(''.join([random.choice("ACGT") for _ in xrange(int(sys.argv[1]))]))
elif len(sys.argv) == 6:
    choices = [('A',float(sys.argv[2])), ('C',float(sys.argv[3])), ('G',float(sys.argv[4])), ('T',float(sys.argv[5]))]
    sys.stdout.write(''.join([weightedChoice(choices) for _ in xrange(int(sys.argv[1]))]))
else:
    print "ERROR: Incorrect number of arguments"
    print "USAGE: python randomDNA.py <length> [P(A) P(C) P(G) P(T)]"
    exit(-1)