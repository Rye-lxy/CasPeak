#! /usr/bin/env python3

def reverseComplement(seq):
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]

def sequenceTrimmer(fastaReader, seqALignmentDict, insertFullSeq, maxTrimmerLength, targetStart, targetEnd, padding=0):
    # trim the reads based on the alignment to the insertion sequence
    for seqName, seq in fastaReader:
        if seqName not in seqALignmentDict:
            continue
        targetAln = seqALignmentDict[seqName][0]
        if targetAln.queryStart > maxTrimmerLength:
            continue

        targetStart = targetStart - padding
        targetEnd = targetEnd + padding
        mid = (targetStart + targetEnd) // 2

        if targetAln.queryStrand == "+":
            if not targetStart <= targetAln.refStart <= targetEnd:
                continue
            else:
                seq = insertFullSeq[mid : targetAln.refStart] + seq[targetAln.queryStart:]
                yield seqName, (seq, "+")
        else:
            if not targetStart <= targetAln.refEnd <= targetEnd:
                continue
            else:
                seq = reverseComplement(seq[targetAln.queryStart:]) + insertFullSeq[targetAln.refEnd : mid]
                yield seqName, (seq, "-")

        




        