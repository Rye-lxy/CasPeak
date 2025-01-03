#! /usr/bin/env python3

def peakDetect(bedgraph, minCov, genome):
    peakChr = None
    peakStart = None
    peakEnd = None
    peakCov = 0
    for line in bedgraph:
        fields = line.split()
        if not fields:
            continue
        chr = fields[0]
        start = int(fields[1])
        end = int(fields[2])
        cov = int(fields[3])

        if cov == 0:
            if peakChr is not None and peakCov >= minCov:
                mid = (peakStart + peakEnd) // 2
                yield peakChr, max(0, mid - 200), min(genome[peakChr], mid + 200), peakCov
            peakChr = None
            peakStart = None
            peakEnd = None
            peakCov = 0
            continue

        if peakChr is None: 
            peakChr = chr
            peakStart = start
            peakEnd = end
            peakCov = cov
        elif peakChr != chr:
            if peakCov >= minCov:
                mid = (peakStart + peakEnd) // 2
                yield peakChr, max(0, mid - 200), min(genome[peakChr], mid + 200), peakCov
            peakChr = chr
            peakStart = start
            peakEnd = end
            peakCov = cov
        else:
            if cov > peakCov:
                peakStart = start
                peakEnd = end
                peakCov = cov
    
    if peakChr is not None and peakCov >= minCov:
        mid = (peakStart + peakEnd) // 2
        yield peakChr, max(0, mid - 200), min(genome[peakChr], mid + 200), peakCov


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Peak detector for bedgraph file")
    parser.add_argument("bedgraph", help="Input bedgraph file")
    parser.add_argument("--minCov", type=int, help="Minimum coverage", default=1)
    args = parser.parse_args()

    with open(args.bedgraph, "r") as bedgraph:
        for peakChr, peakStart, peakEnd, peakCov in peakDetect(bedgraph, args.minCov):
            # bed format
            print(peakChr, peakStart, peakEnd, "{}:{}-{}".format(peakChr, peakStart, peakEnd), peakCov, sep="\t")
