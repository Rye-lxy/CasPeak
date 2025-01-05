#! /usr/bin/env python3

from datetime import datetime

from trimmer import reverseComplement

def vcfHeader():
    return f"""##fileformat=VCFv4.2
##fileDate={datetime.now().strftime("%Y%m%d")}
##source=Caspeak
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of structural variant">
##INFO=<ID=matching_classes,Number=.,Type=String,Description="Repeat class names">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE"""

def vcfRecord(chrom, refpos, altStart, altEnd, strand, seq, id):
    alt = seq[altStart-1:altEnd] if strand == "+" else reverseComplement(seq[altStart-1:altEnd])
    ref = alt[0]
    return f"{chrom}\t{refpos}\tcaspeak_{id}\t{ref}\t{alt}\t.\tPASS\tEND={refpos};SVTYPE=INS;SVLEN={altEnd-altStart};matching_classes=LINE/L1\tGT\t./."