#! /usr/bin/env python3

import itertools
from operator import attrgetter
import os
import subprocess

from filter import joinAll
from fileReader import *

def preAssembler(up, down):
    while up or down:
        if up and down:
            upSeq = up.pop()
            downSeq = down.pop()
            title = f">{upSeq[0]}+{downSeq[0]}"
            seq = upSeq[1] + downSeq[1]
        elif up:
            upSeq = up.pop()
            title = f">{upSeq[0]}"
            seq = upSeq[1]
        else:
            downSeq = down.pop()
            title = f">{downSeq[0]}"
            seq = downSeq[1]
        yield f"{title}\n{seq}\n"

def finalAlignmentCheck(alignments, peakChr, peakStart, peakEnd):
    alns = sorted(list(alignments), key=attrgetter("queryStart"))
    if len(alns) < 2:
        return False

    minInsertLen = 300
    shift = 30
    upstreamAln = None
    downstreamAln = None
    insertLen = 0
    checkRange = range(peakStart, peakEnd)
    for aln in alns:
        if aln.refName == peakChr and aln.refEnd in checkRange:
            if downstreamAln is None:
                upstreamAln = aln
                checkRange = range(aln.refEnd-shift, aln.refEnd+shift)
            elif aln.queryStrand == downstreamAln.queryStrand:
                upstreamAln = aln
        if aln.refName == peakChr and aln.refStart in checkRange:
            if upstreamAln is None:
                downstreamAln = aln
                checkRange = range(aln.refStart-shift, aln.refStart+shift)
            elif aln.queryStrand == upstreamAln.queryStrand:
                downstreamAln = aln

        if upstreamAln and downstreamAln:
            if insertLen >= minInsertLen:
                return True
            else:
                upstreamAln = None
                downstreamAln = None
                insertLen = 0
                checkRange = range(peakStart, peakEnd+1)
        elif upstreamAln or downstreamAln:
            if aln is not upstreamAln and aln is not downstreamAln:
                insertLen += aln.getRefLength()
    return False

def validate(*params):
    # params: only 1.(args) or 2.(args, trimmedReads, peaks)
    # 1. args: trim_read, peak_bed, thread, exog
    # 2. args: thread, exog
    args = params[0]
    if len(params) == 3:
        trimmedReads, peaks = params[1], params[2]
        peaks = peaks.split("\n")
    else:
        trimmedReads = dict(fastaReader(openFile(args.trim_read)))
        trimmedReads = {name[:-1]: (trimmedReads[name], name[-1]) for name in trimmedReads}
        peaks = openFile(args.peak_bed)
    
    os.makedirs("result", exist_ok=True)
    resMaf = open("result/validate.maf", "w")
    resPeak = open("result/validate.peak.bed", "w")
    print("# caspeak validated", file=resMaf, end="\n\n")

    count = 1
    for peak in peaks:
        peakChr, peakStart, peakEnd, _, peakCov = peak.split()
        if int(peakCov) > 1000: # skip peaks with extremely high coverage
            continue

        peakBedData = subprocess.run(["bedtools", "intersect", "-wa", "-a", "peak/sorted.bed", "-b", "-"], capture_output=True, check=True,
                                     input=peak.encode()).stdout.decode().rstrip().split("\n")
        seqNames = set(x.split("\t")[3] for x in peakBedData)

        upstreamReads = [(name, trimmedReads[name][0]) for name in seqNames if trimmedReads[name][1] == "-"]
        downstreamReads = [(name, trimmedReads[name][0]) for name in seqNames if trimmedReads[name][1] == "+"]

        if not upstreamReads or not downstreamReads:
            continue

        assembleProc = subprocess.Popen(["lamassemble", "promethion-2019", "-n", f"peak{count}", "-"], 
                                        stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        sedProc = subprocess.Popen(["sed", "/>/!y/acgt/ACGT/"], 
                                   stdin=assembleProc.stdout, stdout=subprocess.PIPE)
        assembleProc.stdout.close()
        alignProc = subprocess.Popen(["lastal", f"-P{args.thread}", "--split", "lastdb/validate", "-"], 
                                     stdin=sedProc.stdout, stdout=subprocess.PIPE)
        sedProc.stdout.close()

        for seq in preAssembler(upstreamReads, downstreamReads):
            assembleProc.stdin.write(seq.encode())
        assembleProc.stdin.close()
        alignedPeakMaf, _ = alignProc.communicate()

        alignedPeakMafContent = alignedPeakMaf.decode().split("\n")
        if not args.exog and not finalAlignmentCheck(mafReader(alignedPeakMafContent), peakChr, int(peakStart), int(peakEnd)):
            continue

        alignedPeakMafContent = [x for x in alignedPeakMafContent if not x.startswith("#")]
        print("\n".join(alignedPeakMafContent), file=resMaf, end="\n")
        print(peak, file=resPeak, end="\n")
        count += 1

    resMaf.close()