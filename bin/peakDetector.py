#! /usr/bin/env python3

def peakDetect(bedgraph, minCov, minWidth):
    peakChr = None
    peakStart = None
    peakEnd = None
    peakCov = None
    for line in bedgraph:
        fields = line.split()
        if not fields:
            continue
        chr = fields[0]
        start = int(fields[1])
        end = int(fields[2])
        cov = int(fields[3])

        if cov < minCov:
            continue
        
        if peakChr is None: 
            peakChr = chr
            peakStart = start
            peakEnd = end
            peakCov = cov
        elif peakChr == chr and start - peakEnd < 200:
            peakEnd = end
            peakCov = max(peakCov, cov)
        else:
            if peakEnd - peakStart >= minWidth:
                yield peakChr, peakStart, peakEnd, peakCov
            peakChr = chr
            peakStart = start
            peakEnd = end
            peakCov = cov
    
    if peakChr is not None and peakEnd - peakStart >= minWidth:
        yield peakChr, peakStart, peakEnd, peakCov


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Peak detector for bedgraph file")
    parser.add_argument("bedgraph", help="Input bedgraph file")
    parser.add_argument("minWidth", type=int, help="Minimum width of a peak")
    args = parser.parse_args()

    with open(args.bedgraph, "r") as bedgraph:
        for peakChr, peakStart, peakEnd, peakCov in peakDetect(bedgraph, 10, args.minWidth):
            # bed format
            print(peakChr, peakStart, peakEnd, "{}:{}-{}".format(peakChr, peakStart, peakEnd), peakCov, sep="\t")
