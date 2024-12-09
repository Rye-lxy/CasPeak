#! /usr/bin/env python3

from operator import attrgetter
import os
import subprocess

from fileReader import *

def preAssembler(up, down, sample):
    count = 1
    while up or down and count <= sample:
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
        count += 1

def operlapLength(intervalStart1, intervalEnd1, intervalStart2, intervalEnd2):
    if intervalStart1 >= intervalEnd2 or intervalEnd1 <= intervalStart2:
        return 0
    return min(intervalEnd1, intervalEnd2) - max(intervalStart1, intervalStart2)

def finalAlignmentCheck(refAlns, insAlns, peakChr, peakStart, peakEnd):
    alns = sorted(list(refAlns), key=attrgetter("queryStart"))
    # insertAlns = sorted(list(insAlns), key=attrgetter("queryStart"))
    if len(alns) < 2:
        return False

    minInsertLen = 300
    shift = 30
    upstreamAln = None
    downstreamAln = None
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
            insertQueryStart = min(upstreamAln.queryEnd, downstreamAln.queryEnd)
            insertQueryEnd = max(upstreamAln.queryStart, downstreamAln.queryStart)
            insertLen = 0
            for insertAln in insAlns:
                insertLen += operlapLength(insertQueryStart, insertQueryEnd, insertAln.queryStart, insertAln.queryEnd)
            if insertLen >= minInsertLen:
                return True
            else:
                upstreamAln = None
                downstreamAln = None
                checkRange = range(peakStart, peakEnd+1)
    return False

def validate(*params):
    # params: only 1.(args) or 2.(args, trimmedReads, peaks)
    # 1. args: trim_read, peak_bed, thread, exog, sample
    # 2. args: thread, exog, sample
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
        for seq in preAssembler(upstreamReads, downstreamReads, args.sample):
            assembleProc.stdin.write(seq.encode())
        assembleProc.stdin.close()

        assemblyFasta, _ = sedProc.communicate()

        alignValidMaf = subprocess.run(["lastal", f"-P{args.thread}", "--split", "lastdb/validate", "-"], check=True, capture_output=True,
                                    input=assemblyFasta).stdout.decode().split("\n")

        alignInsertMaf = subprocess.run(["lastal", f"-P{args.thread}", "--split", "lastdb/insert", "-"], check=True, capture_output=True,
                                     input=assemblyFasta).stdout.decode().split("\n")

        print("\n".join(alignValidMaf), file=sys.stderr)
        print("\n".join(alignInsertMaf), file=sys.stderr)

        if not args.exog and not finalAlignmentCheck(mafReader(alignValidMaf), mafReader(alignInsertMaf), peakChr, int(peakStart), int(peakEnd)):
            continue

        alignValidMaf = [x for x in alignValidMaf if not x.startswith("#")]
        print("\n".join(alignValidMaf), file=resMaf, end="\n")
        print(peak, file=resPeak, end="\n")
        count += 1

    resMaf.close()