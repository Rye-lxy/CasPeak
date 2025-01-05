#! /usr/bin/env python3

import itertools
import os
import shutil
import subprocess

from .fileReader import *
from .alignmentFilter import *
from .trimmer import *
from .peakDetector import *


def peakAnalyze(args):
    # args: read, ref, insert, genome_maf, insert_maf, target_start, target_end, exog, thread, bedtools_genome, mask
    # min_read_length, max_prop, min_prop, max_trim_length, padding, min_cov, min_width, 
    os.makedirs("peak", exist_ok=True)

    _, insertSeq = next(fastaReader(openFile(args.insert)))

    insertAlignments = dict((name, sorted(list(alns), key=attrgetter("queryStart"))) 
                            for name, alns in itertools.groupby(mafReader(openFile(args.insert_maf)), key=attrgetter("queryName")))
    
    genomeAlignments = dict(genomeAlignmentFilter(mafReader(openFile(args.genome_maf)),
                                                  minReadLen=args.min_read_length,
                                                  maxProp=args.max_prop,
                                                  minProp=args.min_prop,
                                                  exogenous=args.exog))

    trimmedReads = dict(sequenceTrimmer(fastaReader(openFile(args.read)), 
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
    
    try:
        # run bedtools for sorting can calculating coverage
        sortedBed = subprocess.run(["bedtools", "sort", "-i", "-"], capture_output=True, check=True,
                                input="".join(str(aln) for aln in genomeAlignments.values()).encode()).stdout
        genomeCov = subprocess.run(["bedtools", "genomecov", "-bga", "-i", "-", "-g", args.bedtools_genome], capture_output=True, check=True,
                                input=sortedBed).stdout.decode().rstrip()

        bedFile = open("peak/sorted.bed", "w")
        print(sortedBed.decode().rstrip(), file=bedFile)
        bedFile.close()

        # find peaks and store in bed format for bedtools
        peakBed = [f"{name}\t{start}\t{end}\t{name}:{start}-{end}\t{cov}" 
                for name, start, end, cov in peakDetect(genomeCov.split("\n"), args.min_cov, genomeReader(openFile(args.bedtools_genome)))]
        
        # prepare the validate lastdb here to avoid replicate parameters
        os.makedirs("lastdb", exist_ok=True)
        if args.exog:
            subprocess.check_call(f"lastdb -P{args.thread} -uRY4 lastdb/validate {args.ref} {args.insert}", shell=True)
            peaks = "\n".join(peakBed).rstrip()
        else:
            if args.mask:
                if shutil.which("seg-import") is None:
                    print("seg-import is not found in PATH", file=sys.stderr)
                    exit(1)
                subprocess.check_call(f"seg-import rmsk {args.mask} | seg-mask -c - {args.ref} | lastdb -P{args.thread} -uRY4 -R11 -c lastdb/validate - {args.insert}", shell=True)
            else:
                # copy the reference lastdb to validate
                for file in os.listdir("lastdb"):
                    if file.startswith("ref"):
                        shutil.copy(f"lastdb/{file}", "lastdb/validate"+file[3:])

            os.makedirs("lastal", exist_ok=True)
            target_fasta = f">target\n{insertSeq[args.target_start:args.target_end+1]}"
            subprocess.run(["lastdb", "lastdb/target"], input=target_fasta.encode(), check=True)
            with open("lastal/target_to_ref.bed", "w") as f:
                subprocess.check_call(f"lastal -P{args.thread} lastdb/target {args.ref} | maf-convert bed -s 2", stdout=f, shell=True)
            
            peaks = subprocess.run(["bedtools", "subtract", "-A", "-a", "-", "-b", "lastal/target_to_ref.bed"], capture_output=True, check=True,
                                input="\n".join(peakBed).encode()).stdout.decode().rstrip()

        return trimmedReads, peaks
    except subprocess.CalledProcessError as e:
        print(e.stderr.decode(), file=sys.stderr)
        exit(1)
