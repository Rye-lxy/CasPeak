#! /usr/bin/env python3

import sys
import gzip
from Alignment import Alignment

def openFile(fileName):
    if fileName == "-":
        return sys.stdin
    if fileName.endswith(".gz"):
        return gzip.open(fileName, "rt")  # xxx dubious for Python2
    return open(fileName)

def mafReader(lines):
    isRef = True
    for line in lines:
        fields = line.rstrip().split()
        if line.startswith("s"):
            if isRef:
                rName = fields[1]
                rStart = int(fields[2])
                rEnd = rStart + int(fields[3])
                rLength = int(fields[5])
                isRef = False
            else:
                qName = fields[1]
                qStrand = fields[4]
                qLength = int(fields[5])
                if qStrand == "+":
                    qStart = int(fields[2])
                    qEnd = qStart + int(fields[3])
                else:
                    qEnd = qLength - int(fields[2])
                    qStart = qEnd - int(fields[3])
                yield Alignment(rName, rLength, rStart, rEnd, qName, qLength, qStart, qEnd, qStrand)
                isRef = True
    
def fastaReader(lines):
    name = None
    seq = ""
    for line in lines:
        if line.startswith(">"):
            if name:
                yield name, seq
            name = line[1:].split()[0]
            seq = ""
        else:
            seq += line.rstrip()
    if name:
        yield name, seq

