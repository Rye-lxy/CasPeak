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
            title = ">{}+{}".format(upSeq[0], downSeq[0])
            seq = upSeq[1] + downSeq[1]
            yield "{}\n{}\n".format(title, seq)
        elif up:
            upSeq = up.pop()
            title = ">{}".format(upSeq[0])
            yield "{}\n{}\n".format(title, upSeq[1])
        else:
            downSeq = down.pop()
            title = ">{}".format(downSeq[0])
            yield "{}\n{}\n".format(title, downSeq[1])

# TODO: need modification
def finalAlignmentCheck(alignments, peakChr, peakStart, peakEnd):
    alns = sorted(list(alignments), key = attrgetter("queryStrand", "refName", "refStart"))
    if len(alns) < 2 or len(alns) > 25:
        return False

    joinedAlns = joinAll(alns, 200, 7000)
    if len(joinedAlns) <= 2:
        return False

    findLoc = False
    peakRange = range(int(peakStart), int(peakEnd))
    for info, alnChunk in itertools.groupby(alns, key = attrgetter("queryStrand", "refName")):
        _, refName = info
        if refName != peakChr:
            continue
        lastEnd = None
        for aln in alnChunk:
            if lastEnd is None:
                lastEnd = aln.refEnd
            elif lastEnd in peakRange and aln.refStart in peakRange:
                findLoc = True
                break
            else:
                lastEnd = aln.refEnd
        if findLoc:
            break
    
    return findLoc

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
    print("# caspeak validated", file=resMaf, end="\n\n")

    count = 1
    for peak in peaks:
        print(peak)
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
        if not args.exog and not finalAlignmentCheck(mafReader(alignedPeakMafContent), peakChr, peakStart, peakEnd):
            continue

        print("")
        alignedPeakMafContent = [x for x in alignedPeakMafContent if not x.startswith("#")]
        print("\n".join(alignedPeakMafContent), file=resMaf, end="\n")
        count += 1

    resMaf.close()