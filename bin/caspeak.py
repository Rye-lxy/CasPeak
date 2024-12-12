#! /usr/bin/env python3

import argparse
import os
import shutil
import sys

from aligner import aligner
from peakAnalyzer import peakAnalyze
from validator import validate
from plotter import lastDotplot

def addSharedArgs(parser: argparse.ArgumentParser):
    parser.add_argument("--workdir", default=os.getcwd(), metavar="DIR", help="working directory (default: current directory)")
    parser.add_argument("--thread", type=int, default=1, help="number of threads (default: %(default)s)")

def addfilterGroup(parser: argparse.ArgumentParser):
    filterGroup = parser.add_argument_group("parameters for filtering reads")
    filterGroup.add_argument("--min-read-length", type=int, default=500, metavar="NUM",
                             help="minimum read length (default: %(default)s)")
    filterGroup.add_argument("--max-prop", type=float, default=0.99, metavar="NUM",
                             help="maximum proportion of the read aligned to the reference genome (default: %(default)s)")
    filterGroup.add_argument("--min-prop", type=float, default=0.4, metavar="NUM",
                            help="minimum proportion of the read aligned to the reference genome (default: %(default)s)")
    
def addreadsGroup(parser: argparse.ArgumentParser):
    readsGroup = parser.add_argument_group("parameters for trimming reads")
    readsGroup.add_argument("--max-trim-length", type=int, default=100, metavar="NUM",
                            help="maximum length of the adaptor to be trimmed (default: %(default)s)")
    readsGroup.add_argument("--padding", type=int, default=20, metavar="NUM",
                            help="padding for the target region (default: %(default)s)")
    
def addpeakGroup(parser: argparse.ArgumentParser):
    peakGroup = parser.add_argument_group("parameters for peak detection")
    peakGroup.add_argument("--min-cov", type=int, default=10, metavar="NUM",
                           help="minimum coverage for a peak (default: %(default)s)")
    peakGroup.add_argument("--min-width", type=int, default=300,metavar="NUM", 
                           help="minimum width for a peak (default: %(default)s)")
    
def addValidGroup(parser: argparse.ArgumentParser):
    validGroup = parser.add_argument_group("parameters for validation")
    validGroup.add_argument("--sample", type=int, default=500, metavar="N",
                            help="at most N reads are selected for assembly (default: %(default)s)")
    validGroup.add_argument("--vcf", action="store_true", help="output the validated result in VCF format")

def softwareCheck():
    if shutil.which("bedtools") is None:
        print("bedtools is not found in PATH", file=sys.stderr)
        exit(1)
    if shutil.which("lamassemble") is None:
        print("lamassemble is not found in PATH", file=sys.stderr)
        exit(1)
    if shutil.which("lastal") is None:
        print("last is not found in PATH", file=sys.stderr)
        exit(1)

def checkWorkdir(args):   
    os.makedirs(args.workdir, exist_ok=True)
    os.chdir(args.workdir)

def execAlign(args):
    softwareCheck()
    checkWorkdir(args)
    aligner(args)

def execPeak(args):
    softwareCheck()
    if not args.bedtools_genome:
        args.bedtools_genome = "/".join(shutil.which("bedtools").split("/")[:-2] + ["genomes", "human.hg38.genome"])
    checkWorkdir(args)

    trimmedReads, peaks = peakAnalyze(args)

    trimmedReadsFile = open("peak/trimmed_reads.fasta", "w")
    for name, seqInfo in trimmedReads.items():
        seq, strand = seqInfo
        print(f">{name+strand}\n{seq}", file=trimmedReadsFile, end="\n")
    trimmedReadsFile.close()

    peaksFile = open("peak/peaks.bed", "w")
    print(peaks, file=peaksFile)
    peaksFile.close()

def execValid(args):
    softwareCheck()
    checkWorkdir(args)
    validate(args)

def execAllSteps(args):
    softwareCheck()
    if not args.bedtools_genome:
        args.bedtools_genome = "/".join(shutil.which("bedtools").split("/")[:-2] + ["genomes", "human.hg38.genome"])
    checkWorkdir(args)

    aligner(args)

    args.genome_maf = "lastal/read_to_ref.maf"
    args.insert_maf = "lastal/read_to_insert.maf"
    trimmedReads, peaks = peakAnalyze(args)
    trimmedReadsFile = open("peak/trimmed_reads.fasta", "w")
    for name, seqInfo in trimmedReads.items():
        seq, strand = seqInfo
        print(f">{name+strand}\n{seq}", file=trimmedReadsFile, end="\n")
    trimmedReadsFile.close()

    peaksFile = open("peak/peaks.bed", "w")
    print(peaks, file=peaksFile)
    peaksFile.close()

    validate(args, trimmedReads, peaks)
    
def execPlot(args):
    args.background_color = "white"
    args.label_space = 5
    dirPath = os.path.dirname(args.prefix)
    os.makedirs(dirPath, exist_ok=True)
    prefix = os.path.basename(args.prefix)
    if not prefix:
        print("prefix is not specified", file=sys.stderr)
        exit(1)
    lastDotplot(args)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Caspeak: outer-Cas9-targeted nanopore sequencing analysis")
    subparsers = parser.add_subparsers(title="subcommands")

    # align subcommand
    align_parser = subparsers.add_parser("align", help="Align reads to the reference genome and consensus sequence")
    align_parser.add_argument("--read", required=True, metavar="FILE", help="the read FASTA/FATSQ file (required)")
    align_parser.add_argument("--ref", required=True, metavar="FILE", help="the reference genome FASTA file (required)")
    align_parser.add_argument("--insert", required=True, metavar="FILE", help="insert consensus sequence FASTA file (required)")
    addSharedArgs(align_parser)
    align_parser.set_defaults(func=execAlign)

    # peak subcommand
    peak_parser = subparsers.add_parser("peak", help="Detect coverage peaks from the alignment files")
    peak_parser.add_argument("--read", required=True, metavar="FILE", help="the read FASTA/FATSQ file (required)")
    peak_parser.add_argument("--ref", required=True, metavar="FILE", help="the reference genome FASTA file (required)")
    peak_parser.add_argument("--insert", required=True, metavar="FILE", help="insert consensus sequence FASTA file (required)")
    peak_parser.add_argument("--genome-maf", default="lastal/read_to_ref.maf", metavar="MAF", help="alignment to the genome (default: %(default)s)")
    peak_parser.add_argument("--insert-maf", default="lastal/read_to_insert.maf", metavar="MAF", help="alignment to the insertion sequence (default: %(default)s)")
    peak_parser.add_argument("--bedtools-genome", metavar="GENOME", help="genome data for bedtools (default: hg38 file in bedtools)")
    peak_parser.add_argument("-x", "--exog", action="store_true", help="set the insertion exogenous")
    peak_parser.add_argument("--target-start", type=int, required=True, metavar="START",
                            help="start position of the target in the insert sequence (required)")
    peak_parser.add_argument("--target-end", type=int, required=True, metavar="END",
                            help="end position of the target in the insert sequence (required)")
    peak_parser.add_argument("--mask", metavar="RMSK", help="rmsk file on the reference genome generated by RepeatMasker")
    addSharedArgs(peak_parser)
    addfilterGroup(peak_parser)
    addreadsGroup(peak_parser)
    addpeakGroup(peak_parser)
    peak_parser.set_defaults(func=execPeak)

    # valid subcommand
    valid_parser = subparsers.add_parser("valid", help="Validate peaks and output the final result")
    valid_parser.add_argument("--trim-read", required=True, metavar="FILE", help="read FASTA file after trimming (required)")
    valid_parser.add_argument("--peak-bed", required=True, metavar="FILE", help="peak BED file (required)")
    addSharedArgs(valid_parser)
    addValidGroup(valid_parser)
    valid_parser.set_defaults(func=execValid)

    # exec subcommand
    exec_parser = subparsers.add_parser("exec", help="Execute align, peak and valid steps of caspeak")
    exec_parser.add_argument("--read", required=True, metavar="FILE", help="the read FASTA/FATSQ file (required)")
    exec_parser.add_argument("--ref", required=True, metavar="FILE", help="the reference genome FASTA file (required)")
    exec_parser.add_argument("--insert", required=True, metavar="FILE", help="insert consensus sequence FASTA file (required)")
    exec_parser.add_argument("--target-start", type=int, required=True, metavar="START",
                            help="start position of the target in the insert sequence (required)")
    exec_parser.add_argument("--target-end", type=int, required=True, metavar="END",
                            help="end position of the target in the insert sequence (required)")
    exec_parser.add_argument("-x", "--exog", action="store_true", help="set the insertion exogenous")
    exec_parser.add_argument("--bedtools-genome", metavar="GENOME", help="genome data for bedtools (default: hg38 file in bedtools)")
    exec_parser.add_argument("--mask", metavar="RMSK", help="rmsk file on the reference genome generated by RepeatMasker")
    addSharedArgs(exec_parser)
    addfilterGroup(exec_parser)
    addreadsGroup(exec_parser)
    addpeakGroup(exec_parser)
    addValidGroup(exec_parser)
    exec_parser.set_defaults(func=execAllSteps)

    # plot subcommand
    plot_parser = subparsers.add_parser("plot", help="Plot the results from maf file, based on last-dotplot")
    plot_parser.add_argument("--maf", required=True, help="the maf file to plot (required)")
    plot_parser.add_argument("--prefix", default="fig/peak", help="prefix of the output figure name (default: %(default)s)")
    plot_parser.add_argument("-v", "--verbose", action="count",
                  help="show progress messages & data about the plot")
    plot_parser.add_argument("-x", "--width", metavar="INT", type=int, default=1000,
                  help="maximum width in pixels (default: %(default)s)")
    plot_parser.add_argument("-y", "--height", metavar="INT", type=int, default=1000,
                  help="maximum height in pixels (default: %(default)s)")
    plot_parser.add_argument("-m", "--maxseqs", type=int, default=100, metavar="M",
                  help="maximum number of horizontal or vertical sequences "
                  "(default=%(default)s)")
    plot_parser.add_argument("-1", "--seq1", metavar="PATTERN", action="append",
                  default=[],
                  help="which sequences to show from the 1st genome")
    plot_parser.add_argument("-2", "--seq2", metavar="PATTERN", action="append",
                  default=[],
                  help="which sequences to show from the 2nd genome")
    plot_parser.add_argument("--alignments", metavar="FILE", help="secondary alignments")
    plot_parser.add_argument("--sort1", default="1", metavar="N",
                  help="genome1 sequence order: 0=input order, 1=name order, "
                  "2=length order, 3=alignment order (default=%(default)s)")
    plot_parser.add_argument("--sort2", default="1", metavar="N",
                  help="genome2 sequence order: 0=input order, 1=name order, "
                  "2=length order, 3=alignment order (default=%(default)s)")
    plot_parser.add_argument("--strands1", default="0", metavar="N", help=
                  "genome1 sequence orientation: 0=forward orientation, "
                  "1=alignment orientation (default=%(default)s)")
    plot_parser.add_argument("--strands2", default="0", metavar="N", help=
                  "genome2 sequence orientation: 0=forward orientation, "
                  "1=alignment orientation (default=%(default)s)")
    plot_parser.add_argument("--max-gap1", metavar="FRAC", default="1,4", help=
                  "maximum unaligned (end,mid) gap in genome1: "
                  "fraction of aligned length (default=%(default)s)")
    plot_parser.add_argument("--max-gap2", metavar="FRAC", default="1,4", help=
                  "maximum unaligned (end,mid) gap in genome2: "
                  "fraction of aligned length (default=%(default)s)")
    plot_parser.add_argument("--pad", metavar="FRAC", type=float, default=0.04, help=
                  "pad length when cutting unaligned gaps: "
                  "fraction of aligned length (default=%(default)s)")
    plot_parser.add_argument("-j", "--join", default="0", metavar="N", help=
                  "join: 0=nothing, 1=alignments adjacent in genome1, "
                  "2=alignments adjacent in genome2 (default=%(default)s)")
    plot_parser.add_argument("--border-pixels", metavar="INT", type=int, default=1,
                  help="number of pixels between sequences (default=%(default)s)")
    plot_parser.add_argument("-a", "--bed1", "--rmsk1", "--genePred1", "--gap1",
                  action="append", default=[], metavar="FILE",
                  help="read genome1 annotations")
    plot_parser.add_argument("-b", "--bed2", "--rmsk2", "--genePred2", "--gap2",
                  action="append", default=[], metavar="FILE",
                  help="read genome2 annotations")

    textGroup = plot_parser.add_argument_group("Text options")
    textGroup.add_argument("-f", "--fontfile", metavar="FILE",
                  help="TrueType or OpenType font file")
    textGroup.add_argument("-s", "--fontsize", metavar="SIZE", type=int, default=14,
                  help="TrueType or OpenType font size (default: %(default)s)")
    textGroup.add_argument("--labels1", type=int, default=0, metavar="N", help=
                  "genome1 labels: 0=name, 1=name:length, "
                  "2=name:start:length, 3=name:start-end (default=%(default)s)")
    textGroup.add_argument("--labels2", type=int, default=0, metavar="N", help=
                  "genome2 labels: 0=name, 1=name:length, "
                  "2=name:start:length, 3=name:start-end (default=%(default)s)")
    textGroup.add_argument("--rot1", metavar="ROT", default="h",
                  help="text rotation for the 1st genome (default=%(default)s)")
    textGroup.add_argument("--rot2", metavar="ROT", default="v",
                  help="text rotation for the 2nd genome (default=%(default)s)")

    colorGroup = plot_parser.add_argument_group("Color options")
    colorGroup.add_argument("-c", "--forwardcolor", metavar="COLOR", default="red",
                  help="color for forward alignments (default: %(default)s)")
    colorGroup.add_argument("-r", "--reversecolor", metavar="COLOR", default="blue",
                  help="color for reverse alignments (default: %(default)s)")
    colorGroup.add_argument("--border-color", metavar="COLOR", default="black",
                  help="color for pixels between sequences (default=%(default)s)")
    # --break-color and/or --break-pixels for intra-sequence breaks?
    colorGroup.add_argument("--margin-color", metavar="COLOR", default="#dcdcdc",
                  help="margin color")
    colorGroup.add_argument("--exon-color", metavar="COLOR", default="PaleGreen",
                  help="color for exons (default=%(default)s)")
    colorGroup.add_argument("--cds-color", metavar="COLOR", default="LimeGreen",
                  help="color for protein-coding regions (default=%(default)s)")
    colorGroup.add_argument("--bridged-color", metavar="COLOR", default="yellow",
                  help="color for bridged gaps (default: %(default)s)")
    colorGroup.add_argument("--unbridged-color", metavar="COLOR", default="orange",
                  help="color for unbridged gaps (default: %(default)s)")
    
    plot_parser.set_defaults(func=execPlot)

    args = parser.parse_args()

    if hasattr(args, "func"):
        args.func(args)
    else:
        parser.print_help()
        exit(1)