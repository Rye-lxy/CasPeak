#! /usr/bin/env python3

import itertools
from operator import attrgetter

def joinAll(alignments, refGap, queryGap):
    # pick the most representative alignment for one read
    # that indicate the location of insertion
    alns = sorted(list(alignments), key = attrgetter("queryStrand", "refName", "refStart"))
    joinedAlns = []
    for aln in alns:
        if not joinedAlns:
            joinedAlns.append(aln)
            continue

        prevAln = joinedAlns[-1]
        tmpAln = prevAln.join(aln, refGap, queryGap)
        if tmpAln:
            joinedAlns[-1] = tmpAln
        else:
            joinedAlns.append(aln)
    return joinedAlns

def genomeAlignmentFilter(genomeMafReader, minReadLen, maxProp, minProp, exogenous, genomeList):
    # filter for the read alignment to hg38
    for read, alns in itertools.groupby(genomeMafReader, key=attrgetter("queryName", "queryLength")):
        name, length = read
        if length < minReadLen:
            continue
        
        joinedAlns = joinAll(alns, 200, 1000)
        joinedAlns = [aln for aln in joinedAlns if aln.refName in genomeList]
        if not exogenous and len(joinedAlns) < 2:
            continue

        maxAln = max(joinedAlns, key = lambda x: x.getRefLength())
        maxAlnLength = maxAln.getQueryLength()
        if maxAlnLength / length >= maxProp or maxAlnLength / length <= minProp:
            continue
        if length - maxAlnLength < 100:
            continue
        
        maxAln.shrink(toLength=200)
        yield name, maxAln


    