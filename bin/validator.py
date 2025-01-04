#! /usr/bin/env python3

from operator import attrgetter
import os
import subprocess

from fileReader import *
from vcfFormatter import *

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

def overlapLength(intervalStart1, intervalEnd1, intervalStart2, intervalEnd2):
    if intervalStart1 >= intervalEnd2 or intervalEnd1 <= intervalStart2:
        return 0
    return min(intervalEnd1, intervalEnd2) - max(intervalStart1, intervalStart2)

def finalAlignmentCheck(refAlns, insAlns, peakChr, peakStart, peakEnd, minInsertProp):
    alns = sorted(list(refAlns), key=attrgetter("queryStart"))
    insertAlns = list(insAlns)
    if len(alns) < 2:
        return None

    upstreamAln = None
    downstreamAln = None
    pairedAlnsFwd = [None, None]
    pairedAlnsRev = [None, None]
    checkRange = range(peakStart, peakEnd)
    for aln in alns:
        if aln.refName == peakChr and aln.refEnd in checkRange:
            if aln.queryStrand == "+":
                pairedAlnsFwd[0] = aln
            else:
                pairedAlnsRev[0] = aln
        if aln.refName == peakChr and aln.refStart in checkRange:
            if aln.queryStrand == "+":
                pairedAlnsFwd[1] = aln
            else:
                pairedAlnsRev[1] = aln

        if None not in pairedAlnsFwd:
            upstreamAln = pairedAlnsFwd[0]
            downstreamAln = pairedAlnsFwd[1]

            insertQueryStart = min(upstreamAln.queryEnd, downstreamAln.queryEnd)
            insertQueryEnd = max(upstreamAln.queryStart, downstreamAln.queryStart)
            insertLen = 0
            for insertAln in insertAlns:
                insertLen += overlapLength(insertQueryStart, insertQueryEnd, insertAln.queryStart, insertAln.queryEnd)
            if insertLen > 0 and insertLen / (insertQueryEnd - insertQueryStart) >= minInsertProp:
                return upstreamAln.refEnd, insertQueryStart, insertQueryEnd, upstreamAln.queryStrand
            else:
                pairedAlnsFwd = [None, None]
                if aln.refEnd in checkRange and downstreamAln is aln:
                    pairedAlnsFwd[0] = aln
                elif aln.refStart in checkRange and upstreamAln is aln:
                    pairedAlnsFwd[1] = aln
        
        if None not in pairedAlnsRev:
            upstreamAln = pairedAlnsRev[0]
            downstreamAln = pairedAlnsRev[1]

            insertQueryStart = min(upstreamAln.queryEnd, downstreamAln.queryEnd)
            insertQueryEnd = max(upstreamAln.queryStart, downstreamAln.queryStart)
            insertLen = 0
            for insertAln in insertAlns:
                insertLen += overlapLength(insertQueryStart, insertQueryEnd, insertAln.queryStart, insertAln.queryEnd)
            if insertLen > 0 and insertLen / (insertQueryEnd - insertQueryStart) >= minInsertProp:
                return upstreamAln.refStart, insertQueryStart, insertQueryEnd, upstreamAln.queryStrand
            else:
                pairedAlnsRev = [None, None]
                if aln.refEnd in checkRange and downstreamAln is aln:
                    pairedAlnsRev[0] = aln
                elif aln.refStart in checkRange and upstreamAln is aln:
                    pairedAlnsRev[1] = aln

    return None

def validate(*params):
    # params: only 1.(args) or 2.(args, trimmedReads, peaks)
    # 1. args: trim_read, peak_bed, thread, sample, vcf
    # 2. args: thread, sample, vcf
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
    if args.vcf:
        resVcf = open("result/validate.vcf", "w")
        print(vcfHeader(), file=resVcf, end="\n")

    try:
        count = 1
        last = None
        for peak in peaks:
            peakChr, peakStart, peakEnd, _, _ = peak.split()

            peakBedData = subprocess.run(["bedtools", "intersect", "-wa", "-a", "peak/sorted.bed", "-b", "-"], capture_output=True, check=True,
                                        input=peak.encode()).stdout.decode().rstrip().split("\n")
            seqNames = set(x.split("\t")[3] for x in peakBedData)

            upstreamReads = [(name, trimmedReads[name][0]) for name in seqNames if trimmedReads[name][1] == "-"]
            downstreamReads = [(name, trimmedReads[name][0]) for name in seqNames if trimmedReads[name][1] == "+"]

            if not upstreamReads or not downstreamReads:
                continue

            assembleProc = subprocess.Popen(["lamassemble", "promethion-2019", "-n", f"peak{count}", "-P", str(args.thread), "-"], 
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

            vcfData = finalAlignmentCheck(mafReader(alignValidMaf), mafReader(alignInsertMaf), peakChr, int(peakStart), int(peakEnd), args.min_insert)
            if vcfData is None:
                continue

            if last is None:
                last = (peakChr, vcfData[0])
            elif last[0] == peakChr and abs(vcfData[0] - last[1]) < 50:
                continue
            last = (peakChr, vcfData[0])

            alignValidMaf = [x for x in alignValidMaf if not x.startswith("#")]
            print("\n".join(alignValidMaf), file=resMaf, end="\n")
            
            if args.vcf:
                assemblySeq = next(fastaReader(assemblyFasta.decode().split("\n")))[1]
                print(vcfRecord(peakChr, *vcfData, assemblySeq, count), file=resVcf, end="\n")
            
            count += 1

        resMaf.close()
        if args.vcf: resVcf.close()
    except subprocess.CalledProcessError as e:
        print(e.stderr.decode(), file=sys.stderr)
        exit(1)