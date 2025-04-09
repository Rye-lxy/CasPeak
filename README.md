# CasPeak
This is a pipeline for finding non-reference mobile element insertions (MEIs) based on outer-Cas9 targeted Nanopore sequencing and peak detection.

## Installation
You can install CasPeak (together with all requirements) from [Bioconda](https://bioconda.github.io/):
```
conda install -c bioconda caspeak
```
You can also download source code by `git clone` (latest version) or from [release](https://github.com/Rye-lxy/CasPeak/releases).

## Usage
`caspeak` consists of several subcommands:
* [align](#align)
* [peak](#peak)
* [valid](#valid)
* [exec](#exec)
* [plot](#plot)

You can use `caspeak` like:
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

And of course, you should know the location of the Cas9 target site in the consensus sequence.

If you don't care about the intermediate files and procedures, `exec` subcommand provides a shortcut containing all the functions for you:
```
cd /your/work/dir

caspeak exec \
    --read /path/to/read.fa.gz \
    --ref /path/to/reference.fa \
    --insert /path/to/consensus_seq.fa \
    --target-start <NUM> \
    --target-end <NUM> \
    --thread <NUM>  # set the number of parallel CPUs \
    --vcf  # addtional output in VCF format
```

Also, you can run `align`, `peak`, and `valid` sequentially like this script:
```
cd /your/work/dir

caspeak align \
    --read /path/to/read.fa.gz \
    --ref /path/to/reference.fa \
    --insert /path/to/consensus_seq.fa \
    --thread <<NUM>>

caspeak peak \
    --read /path/to/read.fa.gz \
    --ref /path/to/reference.fa \
    --insert /path/to/consensus_seq.fa \
    --target-start <NUM> \
    --target-end <NUM> \
    --thread <NUM>

caspeak valid --thread <NUM> --vcf
```
After running `valid` or `exec`, the final result is named with prefix *validate* in *peak* dir.

Finally, `plot` subcommand can make a dotplot for each peak according to the result MAF file:
```
caspeak plot --maf result/validate.maf
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
| `-v`, `--verbose` | Print more progress messages and data to stdin. |

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
| `--mask` <FILE\> | The region to be masked in either RepeatMasker .out format or UCSC's rmsk format. You can use it to make caspeak execute faster and control the number of false-positive cases, but the false-negative cases might be a bit more. |
| `--thread` <INT\> | Specify the threads running in parallel (default: 1). |
| `--workdir` <DIR\> | Specify the working directory for caspeak output (default: current directory). It is recommended that all the commands should be executed in the same directory. |
||**Parameters for filtering reads**|
| `--min-read-length` <INT\> | Specify the minimum read length N (default: 500). The reads less than N bases will be filtered out. |
| `--max-prop` <NUM\> | Specify the maximum proportion of the longest continuous alignment on the reference genome to the whole read (default: 0.99). | 
| `--min-prop` <NUM\> |  Specify the minimum proportion of the longest continuous alignment on the reference genome to the whole read (default: 0.4). |
||**Parameters for trimming reads**|
| `--max-trim-length` <INT\> | Specify the maximum length N to be trimmed before the first alignment to the insert sequence (default: 100). Reads with distance between the first base and the first alignment longer than N will be filtered out. |
| `--padding` <INT\> | Specify the padding size of the target region (default: 20). This option enables the reads that are not exactly targeted. |
||**Parameters for peak detection**|
| `--min-cov` <INT\> | Specify the minimum coverage to be considered in peak detection (default: 2). |
| `-v`, `--verbose` | Print more progress messages and data to stdin. |
### valid

`caspeak valid` validates peaks listed in BED format. For each peak, `caspeak valid` collects the reads involved, assembles them by [lamassemble](https://gitlab.com/mcfrith/lamassemble), and re-aligns the assembly sequence to validate whether the peak indicates a real non-reference insertion. And it generates the final result to *result* dir only in MAF format if `--vcf` is ignored.
| Option | Description |
| --- | --- |
| `--trim-read` <FILE\> | The read file after filtering and trimming, is usually generated by `caspeak peak` and named as *peak/trimmed_reads.fasta* under your working directory. You can ignore it if the same `--workdir` is specified in `caspeak peak` and `caspeak valid`. |
| `--peak-bed` <FILE\> | The peak file in BED format, usually generated by `caspeak peak` and named as *peak/peaks.bed* under your working directory. You can ignore it if the same `--workdir` is specified in `caspeak peak` and `caspeak valid`. |
| `--thread` <INT\> | Specify the threads running in parallel (default: 1). |
| `--workdir` <DIR\> | Specify the working directory for caspeak output (default: current directory). It is recommended that all the commands should be executed in the same directory. |
| `--sample` <INT\> | Specify that at most N pairs of reads are assembled (default: 20). |
| `--min-insprop` <INT\>| Specify the minimum proportion of the mobile element in the detected insert sequence (default: 0.2). That is, at least 20% of the insert sequence should derive from the mobile element by default. |
| `--min-inslen` <INT\>| Specify the minimum length of the mobile element in the detected insert sequence (default: 80). That is, the length of mobile element in the insert sequence should be longer than 80 bp by default. |
| `--lib` <LIB\> | Use a sequence set of mobile element ancestral lineage (FASTA/FASTAQ format) for validation, instead of a single mobile element sequence specified by `--insert` in previous steps. For example, a set of sequences containing L1HS, L1PA2, L1PA3, ..., L1MA1, ... is suitable for L1HS detection with `--names L1HS`. Some lib files have been prepared in `lib` directory. |
| `--names` <NAME\>| In the LIB file, only the sequences specified in this option are treated as real mobile element for `--min-insert` calculation. Multiple sequence names can be assigned like `--names A --names B --names C`. |
| `--names-re` <REGEXP\> | Treat all the sequences in `lib` file with name matching [REGEXP](https://en.wikipedia.org/wiki/Regular_expression) as real mobile element. Either `--names` or `--names-re` should be set if `--lib` option is specified. |
| `--vcf` | Indicate an extra output in VCF format. |
| `-v`, `--verbose` | Print more progress messages and data to stdin. |
### exec
`caspeak exec` actually provided a shortcut and wrapper for `caspeak align`, `caspeak peak` and `caspeak valid`. It improves speed by skipping several I/O operations.
<table>
    <tr>
        <th scope="col">Option</th>
        <th scope="col">Description</th>
    </tr>
    <tr>
        <td><code>--read</code></td>
        <td rowspan="16">See <a href="#peak">peak options</a>.</td>
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
        <td><code>--sample</code></td>
        <td rowspan="7">See <a href="#valid">valid options</a>.</td>
    </tr>
    <tr>
        <td><code>--min-insprop</code></td>
    </tr>
    <tr>
        <td><code>--min-inslen</code></td>
    </tr>
    <tr>
        <td><code>--lib</code></td>
    </tr>
    <tr>
        <td><code>--names</code></td>
    </tr>
    <tr>
        <td><code>--names-re</code></td>
    </tr>
    <tr>
        <td><code>--vcf</code></td>
    </tr>
    <tr>
        <td><code>-v</code>, <code>--verbose</code></td>
        <td>Print more progress messages and data to stdin.</td>
    </tr>
</table>

### plot
`caspeak plot` makes a dotplot for each peak from the result file in MAF format. It originates from [last-dotplot](https://gitlab.com/mcfrith/last/-/blob/main/doc/last-dotplot.rst) and inherits all options except the input.
| Option | Description |
| --- | --- |
| `--maf` <FILE\> | The result MAF file used for plotting, usually generated by `caspeak valid` or `caspeak exec` (**required**). |
| `--prefix` <PREFIX\> | Specify the name prefix for plots (default: fig/peak). By default, *fig* dir will be made if it doesn't exist, and the plots will be stored as *peak1.png, peak2.png, ...* in *fig* dir. |
| other options | See [last-dotplot document](https://gitlab.com/mcfrith/last/-/blob/main/doc/last-dotplot.rst).|
