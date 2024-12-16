# Caspeak

This is a pipeline for finding non-reference mobile element insertions (MEIs) based on outer-Cas9 targeted Nanopore sequencing and peak detection.

## Usage
`caspeak` consists of several subcommands listed as follows:
* [align](#align)
* [peak](#peak)
* [valid](#valid)
* [exec](#exec)
* [plot](#plot)

You can use `caspeak` as follows:
```
caspeak [subcommand] [options]
```
Specifically, you can check the help message and the options for each subcommand with
```
caspeak [subcommand] -h
```
### Quick Start
Before running caspeak, you need to prepare the following required input file in FASTA or FASTQ format:
* raw Cas9 targeted sequencing data,
* the reference genome (e.g. GRCh38),
* the consensus sequence of the mobile element.

And of course, you should know the where the Cas9 target site is in the consensus sequence.

Then, you can run `align`, `peak`, and `valid` sequentially like this script:
```
caspeak align \
    --read /path/to/read.fa.gz \
    --ref /path/to/reference.fa \
    --insert /path/to/consensus_seq.fa

caspeak peak \
    --read /path/to/read.fa.gz \
    --ref /path/to/reference.fa \
    --insert /path/to/consensus_seq.fa \
    --target-start <NUM> \
    --target-end <NUM>

caspeak valid \
    --trim-read peak/trimmed_reads.fasta \
    --peak-bed peak/peaks.bed \
```
If you don't care about the intermediate files and procedures, `exec` subcommand provides a shortcut containing all the functions for you:
```
caspeak exec \
    --read /path/to/read.fa.gz \
    --ref /path/to/reference.fa \
    --insert /path/to/consensus_seq.fa \
    --target-start <NUM> \
    --target-end <NUM>
```
After running `valid` or `exec`, the final result is named as *validate.maf* in *peak* dir. Using `--vcf` option in either `valid` or `exec` if you need a VCF format output.

Finally, `plot` subcommand can make a dotplot for each peak according to the result MAF file:
```
caspeak plot --maf peak/validate.maf
```
## Subcommands
### align
`caspeak align` utilizes [LAST](https://gitlab.com/mcfrith/last) to align your reads to the reference genome and the insertion sequence. It outputs several Multiple Alignment Format (MAF) files.

| Option | Description |
| --- | --- |
| `--read` <FILE\> | Outer-Cas9 targeted nanopore sequencing data in FASTA or FASTQ format (**required**). |
| `--ref` <FILE\> | Reference genome of the sequenced species in FASTA or FASTQ format (**required**). |
| `--insert` <FILE\> | Consensus sequence of the mobile element in FASTA or FASTQ format (**required**). |
| `--thread` <INT\> | Specify the threads running in parallel (default: 1). |
| `--workdir` <DIR\>| Specify the working directory for caspeak output (default: current directory). It is recommended that all the commands should be executed in the same directory. |

### peak
`caspeak peak` can filter the reads by alignments to both reference genome and insert consensus sequence, and calculate the coverage peaks. Simultaneously, it prepares several files for `caspeak valid` to avoid duplicate parameters.
| Option | Description |
| --- | --- |
| `--read` <FILE\> | Outer-Cas9 targeted nanopore sequencing data in FASTA or FASTQ format (**required**). |
| `--ref` <FILE\> | Reference genome of the sequenced species in FASTA or FASTQ format (**required**). |
| `--insert` <FILE\> | Consensus sequence of the mobile element in FASTA or FASTQ format (**required**). |
| `--target-start` <INT\> | The start position of the Cas9 target site in the consensus sequence (**required**). |
| `--target-end` <INT\> | The end position of the Cas9 target site in the consensus sequence (**required**). |
| `--genome-maf` <FILE\> | The read alignment to the reference genome in MAF format. In general, `caspeak align` outputs it as *lastal/read_to_ref.maf* under your working directory. You can ignore it if you specify the same `--workdir` for `caspeak align` and `caspeak peak`. |
| `--insert-maf` <FILE\> | The read alignment to the consensus sequence in MAF format. In general, `caspeak align` outputs it as *lastal/read_to_insert.maf* and it can also be ignored like `--genome-maf`. |
| `--bedtools-genome` <FILE\> | A genome file passed to [`bedtools genomecov -g`](https://bedtools.readthedocs.io/en/latest/content/tools/genomecov.html). In general, several genome files should be included in *genomes* dir under your path to `bedtools`, and *human.hg38.genome* will be used if you ignore this option. |
| `-x`, `--exog` | Specify the mobile element as an exogenous element, which means there is no identical sequence in the reference genome. |
| `--mask` <FILE\> | The region to be masked in rmsk format. If the reads come from a repeat element (i.e. LINE1) insertion, it is recommended that an annotation file should be sepcified. |
| `--thread` <INT\> | Specify the threads running in parallel (default: 1). |
| `--workdir` <DIR\> | Specify the working directory for caspeak output (default: current directory). It is recommended that all the commands should be executed in the same directory. |
||**Parameters for filtering reads**|
| `--min-read-length` <INT\> | Specify the minimum read leangth N (default: 500). The reads less than N bases will be filtered out. |
| `--max-prop` <NUM\> | Specify the maximum proportion of the longest continuous alignment on the reference genome to the whole read (default: 0.99). | 
| `--min-prop` <NUM\> |  Specify the minimum proportion of the longest continuous alignment on the reference genome to the whole read (default: 0.4). |
||**Parameters for trimming reads**|
| `--max-trim-length` <INT\> | Specify the maximum length N to be trimmed before the first alignment to in insert sequence (default: 100). Reads with distance between the first base and the the first alignment longer than N will be filtered out. |
| `--padding` <INT\> | Specify the padding size of the target region (default: 20). This option enables the reads that are not exactly targeted. |
||**Parameters for peak detection**|
| `--min-cov` <INT\> | Specify the minimum coverage to be considered in peak detection (default: 10). |
| `--min-width` <INT\> | Specify the minimum width of a peak (default: 300). |
### valid

`caspeak valid` validates peaks listed in BED format. For each peak, `caspeak valid` collects the reads involved, assembles them by [lamassemble](https://gitlab.com/mcfrith/lamassemble), and re-aligns the assembly sequence to validate whether the peak indicates a real non-reference insertion. And it generates the final result to *result* dir only in MAF format if `--vcf` is ignored.
| Option | Description |
| --- | --- |
| `--trim-read` <FILE\> | The read file after filtering and trimming, usually generated by `caspeak peak` and named as *peak/trimmed_reads.fasta* under your working directory (**required**). |
| `--peak-bed` <FILE\> | The peak file in BED format, usually generated by `caspeak peak` and named as *peak/peaks.bed* under your working directory (**required**). |
| `--thread` <INT\> | Specify the threads running in parallel (default: 1). |
| `--workdir` <DIR\> | Specify the working directory for caspeak output (default: current directory). It is recommended that all the commands should be executed in the same directory. |
| `--sample` <INT\> | Specify that at most N pairs of reads are assembled (default: 500). |
| `--vcf` | Indicate an extra output in VCF format. |
### exec
`caspeak exec` actually provided a shortcut and wrapper for `caspeak align`, `caspeak peak` and `caspeak valid`. It improves speed by skipping several I/O operations.
<table>
    <tr>
        <th scope="col">Option</th>
        <th scope="col">Description</th>
    </tr>
    <tr>
        <td><code>--read</code></td>
        <td rowspan="17">See <a href="#peak">peak options</a>.</td>
    </tr>
    <tr>
        <td><code>--ref</code></td>
    </tr>
    <tr>
        <td><code>--insert</code></td>
    </tr>
    <tr>
        <td><code>--target-start</code></td>
    </tr>
    <tr>
        <td><code>--target-end</code></td>
    </tr>
    <tr>
        <td><code>-x</code>, <code>--exog</code></td>
    </tr>
    <tr>
        <td><code>--bedtools-genome</code></td>
    </tr>
    <tr>
        <td><code>--mask</code></td>
    </tr>
    <tr>
        <td><code>--workdir</code></td>
    </tr>
    <tr>
        <td><code>--thread</code></td>
    </tr>
    <tr>
        <td><code>--min-read-length</code></td>
    </tr>
    <tr>
        <td><code>--max-prop</code></td>
    </tr>
    <tr>
        <td><code>--min-prop</code></td>
    </tr>
    <tr>
        <td><code>--max-trim-length</code></td>
    </tr>
    <tr>
        <td><code>--padding</code></td>
    </tr>
    <tr>
        <td><code>--min-cov</code></td>
    </tr>
    <tr>
        <td><code>--min-width</code></td>
    </tr>
    <tr>
        <td><code>--sample</code></td>
        <td rowspan="2">See <a href="#valid">valid options</a>.</td>
    </tr>
    <tr>
        <td><code>--vcf</code></td>
    </tr>
</table>

### plot
`caspeak plot` makes a dotplot, a.k.a. Oxford Grid, for each peak from the result file in MAF format. It originates from [last-dotplot](https://gitlab.com/mcfrith/last/-/blob/main/doc/last-dotplot.rst) and inherits all options except the input.
| Option | Description |
| --- | --- |
| `--maf` <FILE\> | The result MAF file used for plotting, usually generated by `caspeak valid` or `caspeak exec` (**required**). |
| `--prefix` <PREFIX\> | Specify the name prefix for plots (default: fig/peak). By default, *fig* dir will be made if it doesn't exist, and the plots will be stored as *peak1.png, peak2.png, ...* in *fig* dir. |
| other options | See [last-dotplot document](https://gitlab.com/mcfrith/last/-/blob/main/doc/last-dotplot.rst).|
