#! /usr/bin/env python3

import sys
import gzip
from .Alignment import Alignment

def openFile(fileName):
    if fileName is None:
        return []
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
    title = None
    fastq = None
    for line in lines:
        if not line:
            continue
        stripped = line.rstrip()
        if fastq is not None:
            fastq.append(stripped)
            if len(fastq) == 4:
                title = fastq[0][1:].split()[0]
                seq = fastq[1]
                yield title, seq
                fastq = []
        elif line[0] == ">":
            if title:
                yield title, "".join(seq)
            title = stripped[1:].split()[0]
            seq = []
        elif title:
            seq.append(stripped)
        else:
            fastq = [stripped]
    if title:
        yield title, "".join(seq)

def genomeReader(lines):
    genome = {}
    for line in lines:
        fields = line.rstrip().split()
        if not fields:
            continue
        genome[fields[0]] = int(fields[1])
    return genome
