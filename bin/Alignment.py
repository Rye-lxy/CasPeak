#! /usr/bin/env python3

class Alignment:
    def __init__(self, refName, refLength, refStart, refEnd, queryName, queryLength, queryStart, queryEnd, queryStrand):
        self.refName = refName
        self.refLength = refLength
        self.refStart = refStart
        self.refEnd = refEnd
        self.queryName = queryName
        self.queryLength = queryLength
        self.queryStart = queryStart
        self.queryEnd = queryEnd
        self.queryStrand = queryStrand

    def __str__(self):
        # generate bed format string
        return "{}\t{}\t{}\t{}\t{}\t{}\n".format(self.refName, self.refStart, self.refEnd, self.queryName, 0, self.queryStrand)
    
    def getRefLength(self):
        return self.refEnd - self.refStart

    def join(self, other, refGap, queryGap):
        # join two alignments if they are co-linear (probably with small indels)
        # return a new alignment object
        if not isinstance(other, Alignment):
            raise TypeError("An alignment instance must be joined with another alignment instance")
        if self.refName == other.refName and self.queryName == other.queryName and self.queryStrand == other.queryStrand:
            if self.queryStrand == "+":
                if other.queryStart - self.queryEnd <= queryGap and other.refStart - self.refEnd <= refGap:
                    return Alignment(self.refName, self.refLength, self.refStart, other.refEnd, self.queryName, self.queryLength, self.queryStart, other.queryEnd, self.queryStrand)
                elif self.queryStart - other.queryEnd <= queryGap and self.refStart - other.refEnd <= refGap:
                    return Alignment(self.refName, self.refLength, other.refStart, self.refEnd, self.queryName, self.queryLength, other.queryStart, self.queryEnd, self.queryStrand)
            else:
                if other.queryStart - self.queryEnd <= queryGap and self.refStart - other.refEnd <= refGap:
                    return Alignment(self.refName, self.refLength, other.refStart, self.refEnd, self.queryName, self.queryLength, self.queryStart, other.queryEnd, self.queryStrand)
                elif self.queryStart - other.queryEnd <= queryGap and other.refStart - self.refEnd <= refGap:
                    return Alignment(self.refName, self.refLength, self.refStart, other.refEnd, self.queryName, self.queryLength, other.queryStart, self.queryEnd, self.queryStrand)     
        return None
    
    def shrink(self, toLength):
        # shrink the alignment to a certain length
        if self.queryStrand == "+":
            self.queryEnd = self.queryStart + toLength
            self.refEnd = self.refStart + toLength
        else:
            self.queryEnd = self.queryStart + toLength
            self.refStart = self.refEnd - toLength