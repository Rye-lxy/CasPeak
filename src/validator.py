#! /usr/bin/env python3

from operator import attrgetter
from functools import partial
import multiprocessing
import os
import subprocess
import shutil

from .fileReader import *
from .vcfFormatter import *

from .logConfigure import getLogger
logger = getLogger(__name__)

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

def finalAlignmentCheck(refAlns, insertAlns, peakChr, peakStart, peakEnd, minInsertProp, minInsertLen):
    alns = sorted(refAlns, key=attrgetter("queryStart"))
    if len(alns) < 2:
        return None
    insertSeqNames = set(x.refName for x in insertAlns)

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
            insertSeqs = {name: 0 for name in insertSeqNames}
            for insertAln in insertAlns:
                overlapLen = overlapLength(insertQueryStart, insertQueryEnd, insertAln.queryStart, insertAln.queryEnd)
                insertLen += overlapLen
                insertSeqs[insertAln.refName] += overlapLen
            if insertLen > minInsertLen and insertLen / (insertQueryEnd - insertQueryStart) >= minInsertProp:
                return upstreamAln.refEnd, insertQueryStart, insertQueryEnd, upstreamAln.queryStrand, max(insertSeqs, key=insertSeqs.get)
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
            insertSeqs = {name: 0 for name in insertSeqNames}
            for insertAln in insertAlns:
                overlapLen = overlapLength(insertQueryStart, insertQueryEnd, insertAln.queryStart, insertAln.queryEnd)
                insertLen += overlapLen
                insertSeqs[insertAln.refName] += overlapLen
            if insertLen > minInsertLen and insertLen / (insertQueryEnd - insertQueryStart) >= minInsertProp:
                return upstreamAln.refEnd, insertQueryStart, insertQueryEnd, upstreamAln.queryStrand, max(insertSeqs, key=insertSeqs.get)
            else:
                pairedAlnsRev = [None, None]
                if aln.refEnd in checkRange and downstreamAln is aln:
                    pairedAlnsRev[0] = aln
                elif aln.refStart in checkRange and upstreamAln is aln:
                    pairedAlnsRev[1] = aln

    return None

def peakAssemble(args, trimmedReads, peaks):
    os.makedirs("tmp", exist_ok=True)
    with open("tmp/assembly.train", "w") as train:
        assemTrainProc = subprocess.Popen(["last-train", f"-P{args.thread}", "-Q0", "lastdb/ref", "-"], stdin=subprocess.PIPE, stdout=train)
        for name, seqInfo in trimmedReads.items():
            assemTrainProc.stdin.write(f">{name}\n{seqInfo[0]}\n".encode())
        assemTrainProc.stdin.close()
        assemTrainProc.communicate()


    for count, peak in enumerate(peaks, start=1):
        if not peak:
            continue
        peakBedData = subprocess.run(["bedtools", "intersect", "-wa", "-a", "peak/sorted.bed", "-b", "-"], capture_output=True, check=True,
                                    input=peak.encode()).stdout.decode().rstrip().split("\n")
        
        
        seqNames = set(x.split("\t")[3] for x in peakBedData)

        upstreamReads = [(name, trimmedReads[name][0]) for name in seqNames if trimmedReads[name][1] == "-"]
        downstreamReads = [(name, trimmedReads[name][0]) for name in seqNames if trimmedReads[name][1] == "+"]


        logger.debug(f"Peak {count} has {len(upstreamReads)} upstream reads and {len(downstreamReads)} downstream reads")

        if not upstreamReads or not downstreamReads:
            continue
        upstreamReads.sort(key=lambda x: len(x[1]))
        downstreamReads.sort(key=lambda x: len(x[1]))

        validReadNum = min(len(upstreamReads), len(downstreamReads)) * 2
        if validReadNum == 2:
            assemblyFasta = next(preAssembler(upstreamReads, downstreamReads, 1)).split("\n")[1]
            assemblyFasta = f">peak{count}-2\n{assemblyFasta}\n"
            assemblyFasta = assemblyFasta.encode()
        else:
            assembleProc = subprocess.Popen(["lamassemble", "-n", f"peak{count}-{validReadNum}", "-P", str(args.thread), "tmp/assembly.train", "-"], 
                                            stdin=subprocess.PIPE, stdout=subprocess.PIPE)
            sedProc = subprocess.Popen(["sed", "/>/!y/acgt/ACGT/"], 
                                    stdin=assembleProc.stdout, stdout=subprocess.PIPE)
            assembleProc.stdout.close()
            for seq in preAssembler(upstreamReads, downstreamReads, args.sample):
                assembleProc.stdin.write(seq.encode())
            assembleProc.stdin.close()

            assemblyFasta, _ = sedProc.communicate() 
        
        yield peak, validReadNum, assemblyFasta

def validateAssembly(assemblyData, args):
    peak, validReadNum, assemblyFasta = assemblyData
    peakChr, peakStart, peakEnd, _, peakCov = peak.split()
    with open(f"tmp/{peakChr}_{peakStart}_{peakEnd}.train", "w") as train:
        trainReturnCode = subprocess.run(["last-train", "-Q0", "lastdb/validate", "-"], stdout=train, input=assemblyFasta).returncode
    if trainReturnCode != 0:
        return None
    mafData = subprocess.run(["lastal", "-p", f"tmp/{peakChr}_{peakStart}_{peakEnd}.train", "--split", "lastdb/validate", "-"], check=True, capture_output=True,
                                input=assemblyFasta).stdout.decode().split("\n")
    alignValidMaf = list(mafReader(mafData))

    if not args.lib:
        alignInsertMaf = subprocess.run(["lastal", "--split", "lastdb/insert", "-"], check=True, capture_output=True,
                                    input=assemblyFasta).stdout.decode().split("\n")
        alignInsertMaf = list(mafReader(alignInsertMaf))
    else:
        alignInsertMaf = subprocess.run(["lastal", "--split", "lastdb/lib", "-"], check=True, capture_output=True,
                                    input=assemblyFasta).stdout.decode().split("\n")
        alignInsertMaf = list(maf for maf in mafReader(alignInsertMaf) if maf.refName in args.names)

    if not alignValidMaf or not alignInsertMaf:
        return None

    if args.test:
        with open("test.maf", "w") as test:
            print("\n".join(mafData), file=test)

    vcfData = finalAlignmentCheck(alignValidMaf, alignInsertMaf, peakChr, int(peakStart), int(peakEnd), args.min_insprop, args.min_inslen)
    if vcfData is None:
        return None
    
    return peakChr, peakStart, peakEnd, peakCov, validReadNum, assemblyFasta.decode(), mafData, vcfData

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

    try:
        logger.info("Start assembling...")
        assemblyList = list(peakAssemble(args, trimmedReads, peaks))
        
        if args.lib:
            subprocess.check_call(f"lastdb -P{args.thread} -uRY4 lastdb/lib {args.lib}", shell=True)
        
        logger.info(f"Start validating with {args.thread} CPUs...")
        with multiprocessing.Pool(args.thread) as pool:
            resultData = pool.map(partial(validateAssembly, args=args), assemblyList)
        resultData = [x for x in resultData if x is not None]

        os.makedirs("result", exist_ok=True)
        resMaf = open("result/validate.maf", "w")
        print("# caspeak validated", file=resMaf, end="\n\n")
        resBed = open("result/validate.bed", "w")
        resFasta = open("result/validate.fasta", "w")
        if args.vcf:
            resVcf = open("result/validate.vcf", "w")
            print(vcfHeader(), file=resVcf, end="\n")
        
        last = None
        count = 1
        logger.info("Writing results...")
        for peakChr, peakStart, peakEnd, peakCov, validReadNum, assemblyFasta, alignValidMaf, vcfData in resultData:
            if last is None:
                last = (peakChr, vcfData[0])
            elif last[0] == peakChr and abs(vcfData[0] - last[1]) < 100:
                continue
            last = (peakChr, vcfData[0])

            alignValidMaf = [x for x in alignValidMaf if not x.startswith("#")]
            print("\n".join(alignValidMaf), file=resMaf, end="\n")
            print(f"{peakChr}\t{peakStart}\t{peakEnd}\tpeak{count}\t{peakCov}", file=resBed, end="\n")
            print(assemblyFasta.rstrip(), file=resFasta, end="\n")

            if args.vcf:
                assemblySeq = next(fastaReader(assemblyFasta.split("\n")))[1]
                print(vcfRecord(peakChr, *vcfData, assemblySeq, count, validReadNum), file=resVcf, end="\n")
            
            count += 1

        resMaf.close()
        resBed.close()
        resFasta.close()
        if args.vcf: resVcf.close()
        shutil.rmtree("tmp")
        logger.info("Finished.")
    except subprocess.CalledProcessError:
        logger.error("Error in peak validation", exc_info=True)
        shutil.rmtree("tmp")
        exit(1)
