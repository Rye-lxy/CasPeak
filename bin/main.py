#! /usr/bin/env python3

import itertools
import os
import subprocess
import sys

from fileReader import *
from filter import *
from trimmer import *
from peakDetector import *


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

def main(args, plotArgs):
    insertName, insertSeq = next(fastaReader(openFile(args.insert_seq)))

    insertAlignments = dict((name, sorted(list(alns), key=attrgetter("queryStart"))) 
                            for name, alns in itertools.groupby(mafReader(openFile(args.insert_maf)), key=attrgetter("queryName")))
    
    genomeAlignments = dict(genomeAlignmentFilter(mafReader(openFile(args.genome_maf)),
                                                  minReadLen=args.min_read_length,
                                                  maxProp=args.max_prop,
                                                  minProp=args.min_prop,
                                                  exogenous=args.exog))

    trimmedReads = dict(sequenceTrimmer(fastaReader(openFile(args.read_fasta)), 
                                        insertAlignments, 
                                        insertSeq, 
                                        maxTrimmerLength=args.max_trim_length, 
                                        targetStart=args.target_start, 
                                        targetEnd=args.target_end+1, 
                                        padding=args.padding))
    
    # intersection of the reads
    sharedReads = genomeAlignments.keys() & trimmedReads.keys()
    genomeAlignments = {name: genomeAlignments[name] for name in sharedReads}
    trimmedReads = {name: trimmedReads[name] for name in sharedReads}
    
    # run bedtools for sorting can calculating coverage
    bedSortProc = subprocess.Popen(["bedtools", "sort", "-i", "-"], 
                                   stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    sortedBed, _ = bedSortProc.communicate(input="".join(str(aln) 
                                                         for aln in genomeAlignments.values()).encode())
    bedCovProc = subprocess.Popen(["bedtools", "genomecov", "-bg", "-i", "-", "-g", args.bedtools_genome], 
                                  stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    genomeCov, _ = bedCovProc.communicate(input=sortedBed)

    bedFile = open("sorted.bed", "w")
    print(sortedBed.decode().rstrip(), file=bedFile)
    bedFile.close()

    genomeCov = genomeCov.decode().rstrip()
    covFile = open("genomeCov.bg", "w")
    print(genomeCov, file=covFile)
    covFile.close()

    # find peaks and store in bed format for bedtools
    peakBed = [f"{name}\t{start}\t{end}\t{name}:{start}-{end}\t{cov}" 
               for name, start, end, cov in peakDetect(genomeCov.split("\n"), args.min_cov, args.min_width)]
    
    # ignore peaks that overlap with the target regions in hg38
    if args.ignore_bed:
        subtractProc = subprocess.Popen(["bedtools", "subtract", "-A", "-a", "-", "-b", args.ignore_bed], 
                                        stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        peaks, _ = subtractProc.communicate(input="\n".join(peakBed).encode())
        peaks = peaks.decode().rstrip()
    else:
        peaks = "\n".join(peakBed)

    peakFile = open("peaks.bed", "w")
    print(peaks, file=peakFile)
    peakFile.close()

    os.makedirs("fig", exist_ok=True)
    resAlnFile = open("peaks.maf", "a")

    peaks = peaks.rstrip().split("\n")
    count = 0
    for peak in peaks:
        count += 1
        peakChr, peakStart, peakEnd, _, peakCov = peak.split()
        if int(peakCov) > 1000: # skip peaks with extremely high coverage
            continue

        interscetProc = subprocess.Popen(["bedtools", "intersect", "-wa", "-a", "sorted.bed", "-b", "-"], 
                                         stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        cutProc = subprocess.Popen(["cut", "-f", "4"], 
                                   stdin=interscetProc.stdout, stdout=subprocess.PIPE)
        interscetProc.stdout.close()
        sortProc = subprocess.Popen(["sort"], 
                                    stdin=cutProc.stdout, stdout=subprocess.PIPE)
        cutProc.stdout.close()
        uniqProc = subprocess.Popen(["uniq"], 
                                    stdin=sortProc.stdout, stdout=subprocess.PIPE)
        sortProc.stdout.close()

        interscetProc.stdin.write(peak.encode())
        interscetProc.stdin.close()

        seqNames, _ = uniqProc.communicate()
        seqNames = seqNames.decode().rstrip().split("\n")
        upstreamReads = [(name, trimmedReads[name][0]) for name in seqNames if trimmedReads[name][1] == "-"]
        downstreamReads = [(name, trimmedReads[name][0]) for name in seqNames if trimmedReads[name][1] == "+"]

        if not upstreamReads or not downstreamReads:
            continue

        assembleProc = subprocess.Popen(["lamassemble", "promethion-2019", "-n", f"peak{count}", "-"], 
                                        stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        sedProc = subprocess.Popen(["sed", "/>/!y/acgt/ACGT/"], 
                                   stdin=assembleProc.stdout, stdout=subprocess.PIPE)
        assembleProc.stdout.close()
        alignProc = subprocess.Popen(["lastal", "-P"+str(args.thread), "--split", args.lastdb, "-"], 
                                     stdin=sedProc.stdout, stdout=subprocess.PIPE)
        sedProc.stdout.close()

        for seq in preAssembler(upstreamReads, downstreamReads):
            assembleProc.stdin.write(seq.encode())
        assembleProc.stdin.close()
        alignedPeakMaf, _ = alignProc.communicate()

        alignedPeakMafContent = alignedPeakMaf.decode().split("\n")
        if not finalAlignmentCheck(mafReader(alignedPeakMafContent), peakChr, peakStart, peakEnd):
            continue

        alignedPeakMafContent = [x for x in alignedPeakMafContent if not x.startswith("#")]
        print("\n".join(alignedPeakMafContent), file=resAlnFile, end="\n\n")

        plotProc = subprocess.Popen(["last-dotplot"] + plotArgs + ["-", f"fig/peak{count}.png"], stdin=subprocess.PIPE)
        plotProc.communicate(input=alignedPeakMaf)
    
    resAlnFile.close()

if __name__ == "__main__":
    import argparse
    import shutil


    parser = argparse.ArgumentParser()
    parser.add_argument("--genome-maf", required=True, metavar="MAF", help="alignment to the genome (required)")
    parser.add_argument("--insert-maf", required=True, metavar="MAF", help="alignment to the insertion sequence (required)")
    parser.add_argument("--read-fasta", required=True, metavar="FASTA", help="read sequences (required)")
    parser.add_argument("--insert-seq", required=True, metavar="FASTA", help="insertion sequence (required)")
    parser.add_argument("--lastdb", required=True, metavar="LASTDB", help="lastdb for validation genome (required)")
    parser.add_argument("--ignore-bed", metavar="BED", help="regions to exclude from peak detection")
    parser.add_argument("--thread", type=int, default=8, help="number of threads (default: 8)")
    parser.add_argument("--bedtools-genome", metavar="GENOME", help="genome data for bedtools")
    parser.add_argument("-x", "--exog", action="store_true", help="specify the insertion exogenous")
    
    filterGroup = parser.add_argument_group("Arguments for filtering reads")
    filterGroup.add_argument("--min-read-length", type=int, default=500, metavar="NUM",
                             help="minimum read length (default: 500)")
    filterGroup.add_argument("--max-prop", type=float, default=0.99, metavar="NUM",
                             help="maximum proportion of the read aligned to the reference genome (default: 0.99)")
    filterGroup.add_argument("--min-prop", type=float, default=0.4, metavar="NUM",
                            help="minimum proportion of the read aligned to the reference genome (default: 0.4)")
    
    readsGroup = parser.add_argument_group("Arguments for trimming reads")
    readsGroup.add_argument("--max-trim-length", type=int, default=100, metavar="NUM",
                            help="maximum length of the adaptor to be trimmed (default: 100)")
    readsGroup.add_argument("--target-start", type=int, required=True, metavar="START",
                            help="start of the target region (required)")
    readsGroup.add_argument("--target-end", type=int, required=True, metavar="END",
                            help="end of the target region (required)")
    readsGroup.add_argument("--padding", type=int, default=20, metavar="NUM",
                            help="padding for the target region (default: 20)")
    
    peakGroup = parser.add_argument_group("Arguments for peak detection")
    peakGroup.add_argument("--min-cov", type=int, default=10, metavar="NUM",
                           help="minimum coverage for a peak (default: 10)")
    peakGroup.add_argument("--min-width", type=int, default=300,metavar="NUM", 
                           help="minimum width for a peak (default: 300)")

    args, plotArgs = parser.parse_known_args()

    if shutil.which("bedtools") is None:
        print("bedtools is not found in PATH", file=sys.stderr)
        exit(1)
    if shutil.which("lamassemble") is None:
        print("lamassemble is not found in PATH", file=sys.stderr)
        exit(1)
    if shutil.which("lastal") is None:
        print("last is not found in PATH", file=sys.stderr)
        exit(1)

    if not args.bedtools_genome:
        args.bedtools_genome = "/".join(shutil.which("bedtools").split("/")[:-2] + ["genomes", "human.hg38.genome"])

    if not os.path.exists(args.genome_maf):
        print(f"{args.genome_maf} does not exist", file=sys.stderr)
        exit(1)
    if not os.path.exists(args.insert_maf):
        print(f"{args.insert_maf} does not exist", file=sys.stderr)
        exit(1)
    if not os.path.exists(args.read_fasta):
        print(f"{args.read_fasta} does not exist", file=sys.stderr)
        exit(1)
    if not os.path.exists(args.insert_seq):
        print(f"{args.insert_seq} does not exist", file=sys.stderr)
        exit(1)
    if args.ignore_bed and not os.path.exists(args.ignore_bed):
        print(f"{args.ignore_bed} does not exist", file=sys.stderr)
        exit(1)
    if not os.path.exists(args.bedtools_genome):
        print(f"{args.bedtools_genome} does not exist", file=sys.stderr)
        exit(1)

    main(args, plotArgs)
