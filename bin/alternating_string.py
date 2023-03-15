#!/usr/bin/env python3

def is_alternating(string) : 

    if len(string) < 3 : 
        return False

    for i in range(len(string) - 2) : 
        if not string[i] == string[i+2] : 
            return False

    if string[0] == string[1] : 
        return False
    

    return True

if __name__ == "__main__" : 

    tests = ["GTGTGT", "GTGGGT", "G", "GGTGTGT", "GT", "GTG"]
    for t in tests : 
        print( is_alternating(t) )