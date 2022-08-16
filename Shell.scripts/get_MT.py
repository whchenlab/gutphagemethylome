#!/usr/bin/python
# -*- coding: UTF-8 -*-
import argparse
import re

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="mcl file", dest="clustr", type=str, required=True)
parser.add_argument("-o", "--output", help="mcl table", dest="out", type=str, required=True)

args = parser.parse_args()
clustr = args.clustr
outfile = args.out

out = open(outfile,"w")
with open(clustr,"r") as fh :
    i = 1
    for line in fh:
        line = line.strip()
        arr = line.split("\t")
        for a in arr:
            out.write("MTase_"+str(i)+",%s\n"%a)
        i+=1