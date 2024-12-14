# Caspeak

This is a pipeline for finding mobile element insertions (MEIs) based on outer-Cas9 targeted Nanopore sequencing and peak detection.

## Usage
`caspeak` consists of several subcommands, like `git` or `bedtools`, listed as follows:
* [align](#align)
* [peak](#peak)
* [valid](#valid)
* [exec](#exec)
* [plot](#plot)

You can use `caspeak` as follows:
```
caspeak [subcommand] [options]
```
Specifically, you can check the options for each subcommand with
```
caspeak [subcommand] -h
```
### align
`caspeak align` utilizes [LAST](https://gitlab.com/mcfrith/last) to align your reads to the reference genome and the insertion sequence. It will output several Multiple Alignment Format (MAF) files.

| Option | Description |
| --- | --- |
| `--read` | Outer-Cas9 targeted nanopore sequencing data in FASTA or FASTQ format (**required**). |
| `--ref` | Reference genome of the sequenced species in FASTA or FASTQ format (**required**). |
| `--insert` | Consensus sequence of the mobile element in FASTA or FASTQ format (**required**). |
| `--thread` | Specify the threads running in parallel (default: 1). |
| `--workdir`| Specify the root directory for caspeak output (default: current directory). |

### peak
`caspeak peak` can filter the reads by alignments to both reference genome and insert consensus sequence, and calculate the coverage peaks. Simultaneously, it will prepare several files for `caspeak valid` to avoid duplicate parameters.
| Option | Description |
| --- | --- |
| `--read` | Outer-Cas9 targeted nanopore sequencing data in FASTA or FASTQ format (**required**). |
| `--ref` | Reference genome of the sequenced species in FASTA or FASTQ format (**required**). |
| `--insert` | Consensus sequence of the mobile element in FASTA or FASTQ format (**required**). |
| `--target-start` | The start position of the Cas9 target site in the consensus sequence (**required**). |
| `--target-end` | The end position of the Cas9 target site in the consensus sequence (**required**). |
| `--genome-maf` | The read alignment to the reference genome. It can be found automatically if you specify the same `--workdir` in `caspeak align` and `caspeak peak`. |
| `--insert-maf` | The read alignment to the consensus sequence, and it can also be found automatically like `--genome-maf`. |
| `--bedtools-genome` | A genome file passed to [`bedtools genomecov -g`](https://bedtools.readthedocs.io/en/latest/content/tools/genomecov.html). In general, several genome files should be included in the *genomes* directory under your path to `bedtools`, and `caspeak peak` will use *human.hg38.genome* in it if you ignore this option. |
| `-x`, `--exog` | Specify the mobile element as an exogenous element, which means there is no identical sequence in the reference genome. |
| `--mask` | The region to be masked in rmsk format. If the reads come from a repeat element (i.e. LINE1) insertion, it is recommended that an annotation file should be sepcified. |
| `--thread` | specify the threads running in parallel (default: 1) |
| `--workdir`| specify the root directory for caspeak output (default: current directory) |
### valid

### exec

### plot