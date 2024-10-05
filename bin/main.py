#! /usr/bin/env python3

import os
import subprocess
import sys

from fileReader import *
from filter import *
from trimmer import *
from peakDetector import *

BEDTOOLS_GENOME = "/home/rye/.local/bedtools/2.31.1/genomes/human.hg38.genome"
LASTDB = "/big/rye/human/last/R11/softMaskedHg38_L1Hs/lastdb"
LINE_ANNO = "/big/rye/human/seg-mask/segfiles/LINE1.rmsk.txt"

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

def main(args):
    insertName, insertSeq = next(fastaReader(openFile(args.insert_seq)))

    insertAlignments = dict((name, list(alns)) for name, alns in itertools.groupby(mafReader(openFile(args.insert_maf)), key=attrgetter("queryName")))
    genomeAlignments = dict(genomeAlignmentFilter(mafReader(openFile(args.genome_maf))))

    trimmedReads = dict(sequenceTrimmer(fastaReader(openFile(args.read_fasta)), insertAlignments, insertSeq, 
                                        maxTrimmerLength=100, targetStart=3658, targetEnd=3689, padding=20))
    
    # intersection of the reads
    sharedReads = genomeAlignments.keys() & trimmedReads.keys()
    genomeAlignments = {name: genomeAlignments[name] for name in sharedReads}
    trimmedReads = {name: trimmedReads[name] for name in sharedReads}
    
    # run bedtools for sorting can calculating coverage
    bedSortProc = subprocess.Popen(["bedtools", "sort", "-i", "-"], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    sortedBed, _ = bedSortProc.communicate(input="".join(str(aln) for aln in genomeAlignments.values()).encode())
    bedCovProc = subprocess.Popen(["bedtools", "genomecov", "-bg", "-i", "-", "-g", BEDTOOLS_GENOME], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    genomeCov, _ = bedCovProc.communicate(input=sortedBed)

    bedFile = open("sorted.bed", "w")
    print(sortedBed.decode(), file=bedFile)
    bedFile.close()

    # find peaks and store in bed format for bedtools
    peaks = [f"{name}\t{start}\t{end}\t{name}:{start}-{end}\t{cov}" for name, start, end, cov in peakDetect(genomeCov.decode().rstrip().split("\n"), 10, 201)]
    # ignore peaks that overlap with the target regions in hg38
    subtractProc = subprocess.Popen(["bedtools", "subtract", "-A", "-a", "-", "-b", args.ignore_bed], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    peaks, _ = subtractProc.communicate(input="\n".join(peaks).encode())

    os.makedirs("fig", exist_ok=True)
    peaks = peaks.decode().rstrip().split("\n")
    count = 0
    for peak in peaks:
        count += 1
        interscetProc = subprocess.Popen(["bedtools", "intersect", "-wa", "-a", "sorted.bed", "-b", "-"], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        cutProc = subprocess.Popen(["cut", "-f", "4"], stdin=interscetProc.stdout, stdout=subprocess.PIPE)
        interscetProc.stdout.close()
        sortProc = subprocess.Popen(["sort"], stdin=cutProc.stdout, stdout=subprocess.PIPE)
        cutProc.stdout.close()
        uniqProc = subprocess.Popen(["uniq"], stdin=sortProc.stdout, stdout=subprocess.PIPE)
        sortProc.stdout.close()

        interscetProc.stdin.write(peak.encode())
        interscetProc.stdin.close()

        seqNames, _ = uniqProc.communicate()
        seqNames = seqNames.decode().rstrip().split("\n")
        upstreamReads = [(name, trimmedReads[name][0]) for name in seqNames if trimmedReads[name][1] == "-"]
        downstreamReads = [(name, trimmedReads[name][0]) for name in seqNames if trimmedReads[name][1] == "+"]

        if not upstreamReads or not downstreamReads:
            continue

        assembleProc = subprocess.Popen(["lamassemble", "promethion-2019", "-n", f"peak{count}", "-"], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        sedProc = subprocess.Popen(["sed", "/>/!y/acgt/ACGT/"], stdin=assembleProc.stdout, stdout=subprocess.PIPE)
        assembleProc.stdout.close()
        alignProc = subprocess.Popen(["lastal", "-P8", "--split", LASTDB, "-"], stdin=sedProc.stdout, stdout=subprocess.PIPE)
        sedProc.stdout.close()
        plotProc = subprocess.Popen(["last-dotplot", "-a", LINE_ANNO, "--labels1=3", "--strands2=1", "-j", "2", 
                                     "--rot1=v", "-", f"fig/peak{count}.png"], stdin=alignProc.stdout)
        alignProc.stdout.close()
        for seq in preAssembler(upstreamReads, downstreamReads):
            assembleProc.stdin.write(seq.encode())
        assembleProc.stdin.close()
        plotProc.communicate()
        
    os.remove("sorted.bed")

if __name__ == "__main__":
    import argparse
    import shutil

    if shutil.which("bedtools") is None:
        print("bedtools is not found in PATH", file=sys.stderr)
        exit(1)
    if shutil.which("lamassemble") is None:
        print("lamassemble is not found in PATH", file=sys.stderr)
        exit(1)
    if shutil.which("lastal") is None:
        print("last is not found in PATH", file=sys.stderr)
        exit(1)

    parser = argparse.ArgumentParser()
    parser.add_argument("--genome-maf", required=True, help="alignment to the genome")
    parser.add_argument("--insert-maf", required=True, help="alignment to the insertion")
    parser.add_argument("--read-fasta", required=True, help="read sequences")
    parser.add_argument("--insert-seq", required=True, help="insertion sequence")
    parser.add_argument("--ignore-bed", required=True, help="regions to ignore")

    args = parser.parse_args()

    main(args)
