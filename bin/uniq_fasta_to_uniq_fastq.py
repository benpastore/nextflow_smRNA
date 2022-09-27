#!/usr/bin/python

import sys

FILE=sys.argv[1]

with open("{}".format(FILE), "r") as f :
    for line in f:
        if line.startswith(">") : 
            info = line.strip().split(":")
            sequence = info[0].replace(">","")
            count = info[1]
            qualscore = "I"*len(sequence)
            print("@{}:{}\n{}\n+\n{}".format(sequence,count,sequence,qualscore))
f.close()
