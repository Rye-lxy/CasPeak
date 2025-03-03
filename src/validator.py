#! /usr/bin/env python3

from operator import attrgetter
import os
import subprocess

from .fileReader import *
from .vcfFormatter import *

def preAssembler(up, down, sample):
    count = 1
    while up and down and count <= sample:
        upSeq = up.pop()
        downSeq = down.pop()
        title = f">{upSeq[0]}+{downSeq[0]}"
        seq = upSeq[1] + downSeq[1]
        yield f"{title}\n{seq}\n"
        count += 1

def overlapLength(intervalStart1, intervalEnd1, intervalStart2, intervalEnd2):
    if intervalStart1 >= intervalEnd2 or intervalEnd1 <= intervalStart2:
        return 0
    return min(intervalEnd1, intervalEnd2) - max(intervalStart1, intervalStart2)

def finalAlignmentCheck(refAlns, insertAlns, peakChr, peakStart, peakEnd, minInsertProp):
    alns = sorted(list(refAlns), key=attrgetter("queryStart"))
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
            if insertLen > 80 and insertLen / (insertQueryEnd - insertQueryStart) >= minInsertProp:
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
                return upstreamAln.refEnd, insertQueryStart, insertQueryEnd, upstreamAln.queryStrand
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
    resBed = open("result/validate.bed", "w")
    resFasta = open("result/validate.fasta", "w")
    if args.vcf:
        resVcf = open("result/validate.vcf", "w")
        print(vcfHeader(), file=resVcf, end="\n")

    os.makedirs("tmp", exist_ok=True)
    try:
        with open("tmp/assembly.train", "w") as train:
            assemTrainProc = subprocess.Popen(["last-train", f"-P{args.thread}", "-Q0", "lastdb/ref", "-"],
                                            stdin=subprocess.PIPE, stdout=train)
            for name, seqInfo in trimmedReads.items():
                assemTrainProc.stdin.write(f">{name}\n{seqInfo[0]}\n".encode())
            assemTrainProc.stdin.close()
            assemTrainProc.communicate()
    except subprocess.CalledProcessError:
        print("Error in peak validation", file=sys.stderr)
        exit(1)

    count = 1
    last = None
    for peak in peaks:
        peakChr, peakStart, peakEnd, _, peakCov = peak.split()
        try:
            peakBedData = subprocess.run(["bedtools", "intersect", "-wa", "-a", "peak/sorted.bed", "-b", "-"], capture_output=True, check=True,
                                        input=peak.encode()).stdout.decode().rstrip().split("\n")
        except subprocess.CalledProcessError:
            print("Error in peak validation", file=sys.stderr)
            exit(1)
        
        seqNames = set(x.split("\t")[3] for x in peakBedData)

        upstreamReads = [(name, trimmedReads[name][0]) for name in seqNames if trimmedReads[name][1] == "-"]
        downstreamReads = [(name, trimmedReads[name][0]) for name in seqNames if trimmedReads[name][1] == "+"]

        if args.test:
            print(f"Peak {count} has {len(upstreamReads)} upstream reads and {len(downstreamReads)} downstream reads", file=sys.stderr)

        if not upstreamReads or not downstreamReads:
            continue
        upstreamReads.sort(key=lambda x: len(x[1]))
        downstreamReads.sort(key=lambda x: len(x[1]))

        if len(upstreamReads) == 1 or len(downstreamReads) == 1:
            assemblyFasta = next(preAssembler(upstreamReads, downstreamReads, 1)).split("\n")[1]
            assemblyFasta = f">peak{count}\n{assemblyFasta}\n".encode()
        else:
            try:
                assembleProc = subprocess.Popen(["lamassemble", "-n", f"peak{count}", "-P", str(args.thread), "tmp/assembly.train", "-"], 
                                                stdin=subprocess.PIPE, stdout=subprocess.PIPE)
                sedProc = subprocess.Popen(["sed", "/>/!y/acgt/ACGT/"], 
                                        stdin=assembleProc.stdout, stdout=subprocess.PIPE)
                assembleProc.stdout.close()
                for seq in preAssembler(upstreamReads, downstreamReads, args.sample):
                    assembleProc.stdin.write(seq.encode())
                assembleProc.stdin.close()

                assemblyFasta, _ = sedProc.communicate()
            except subprocess.CalledProcessError:
                print("Error in peak validation", file=sys.stderr)
                exit(1)

        try:
            with open("tmp/validate.train", "w") as train:
                subprocess.run(["last-train", f"-P{args.thread}", "-Q0", "lastdb/validate", "-"], check=True, input=assemblyFasta, stdout=train, stderr=subprocess.DEVNULL)
        except subprocess.CalledProcessError:
            continue

        try:
            alignValidMaf = subprocess.run(["lastal", f"-P{args.thread}", "-p", "tmp/validate.train", "--split", "lastdb/validate", "-"], check=True, capture_output=True,
                                        input=assemblyFasta).stdout.decode().split("\n")

            if not args.lib:
                alignInsertMaf = subprocess.run(["lastal", f"-P{args.thread}", "--split", "lastdb/insert", "-"], check=True, capture_output=True,
                                            input=assemblyFasta).stdout.decode().split("\n")
                alignInsertMaf = list(mafReader(alignInsertMaf))
            else:
                subprocess.check_call(f"lastdb -P{args.thread} -uRY4 lastdb/lib {args.lib}", shell=True)
                alignInsertMaf = subprocess.run(["lastal", f"-P{args.thread}", "--split", "lastdb/lib", "-"], check=True, capture_output=True,
                                            input=assemblyFasta).stdout.decode().split("\n")
                alignInsertMaf = list(maf for maf in mafReader(alignInsertMaf) if maf.refName in args.names)
        except subprocess.CalledProcessError:
            print("Error in peak validation", file=sys.stderr)
            exit(1)

        if args.test:
            print("\n".join(alignValidMaf), file=sys.stdout)

        vcfData = finalAlignmentCheck(mafReader(alignValidMaf), alignInsertMaf, peakChr, int(peakStart), int(peakEnd), args.min_insert)        
        if vcfData is None:
            continue

        if last is None:
            last = (peakChr, vcfData[0])
        elif last[0] == peakChr and abs(vcfData[0] - last[1]) < 50:
            continue
        last = (peakChr, vcfData[0])

        alignValidMaf = [x for x in alignValidMaf if not x.startswith("#")]
        print("\n".join(alignValidMaf), file=resMaf, end="\n")
        print(f"{peakChr}\t{peakStart}\t{peakEnd}\tpeak{count}\t{peakCov}", file=resBed, end="\n")
        print(assemblyFasta.decode().rstrip(), file=resFasta, end="\n")

        if args.vcf:
            assemblySeq = next(fastaReader(assemblyFasta.decode().split("\n")))[1]
            print(vcfRecord(peakChr, *vcfData, assemblySeq, count), file=resVcf, end="\n")
        
        count += 1

    resMaf.close()
    resBed.close()
    resFasta.close()
    if args.vcf: resVcf.close()