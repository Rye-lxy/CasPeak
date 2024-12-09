#! /usr/bin/env python3

from datetime import datetime

from trimmer import reverseComplement

def vcfHeader():
    return f"""##fileformat=VCFv4.2
##fileDate={datetime.now().strftime("%Y%m%d")}
##source=Caspeak
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"""

def vcfRecord(chrom, refpos, altStart, altEnd, strand, seq, id):
    alt = seq[altStart-1:altEnd] if strand == "+" else reverseComplement(seq[altStart-1:altEnd])
    ref = alt[0]
    return f"{chrom}\t{refpos}\tcaspeak_{id}\t{ref}\t{alt}\t.\tPASS\t."