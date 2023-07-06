#!/usr/bin/python

import sys
import random

FILE=sys.argv[1]

with open("{}".format(FILE), "r") as f :
    for line in f:
        if line.startswith(">") : 
            info = line.strip()
            sequence = info.replace(">","")
            count = info[1]
            qualscore = "I"*len(sequence)
            print("@{}_x{}\n{}\n+\n{}".format(sequence, random.randint(0,1e9),sequence,qualscore))
f.close()