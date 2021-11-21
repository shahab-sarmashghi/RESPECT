# RESPECT
RESPECT is a tool for estimating the repeat spectra and parameters such as
genome length and coverage from low-coverage sequencing reads (genome-skims). 
The preprint of the manuscript where we describe the algorithm and test 
RESPECT on simulated short reads and SRA data:
 - [Sarmashghi, S., Balaban, M., Rachtman, E., Touri, B., Mirarab, S. and 
    Bafna, V., 2021. Estimating repeat spectra and genome length from 
    low-coverage genome skims with RESPECT. bioRxiv.][1]

RESPECT is a command-line tool implemented in python. It runs [Jellyfish][2] 
internally to efficiently compute k-mer histograms of input sequence files.
RESPECT also uses [Gurobi][3] to accurately solve optimization problems.

Installation
------------
You need to have python 3.6, 3.7, or 3.8 installed. It is recommended to create
a new conda environment and use conda to install the following requirements 
(alternatively, you can install each tool from its source and add them to the 
system path).
1. Install [Miniconda][4] (you can skip this if you already have either of 
Miniconda or Anaconda installed).
2. Add bioconda and Gurobi channels by running the following commands in 
the terminal (order matters):
    ```
        conda config --add channels defaults
        conda config --add channels conda-forge
        conda config --add channels bioconda
        conda config --add channels https://conda.anaconda.org/gurobi
    ```
2. Install [Jellyfish][2], [seqtk][5], and [gurobi][6] by running the following 
command
    ```
        conda install jellyfish seqtk gurobi 
    ```
3. Set up the [license][7] to use Gurobi for academic purposes.
4. (Optional) have `gzip` 1.6 or later installed to process gzipped inputs.
Otherwise, if you have .gz inputs, you need to select a python library 
for decompression using (`--decomp`) option.
5. Clone the github repository by running (or you can download the repo)
    ```
        git clone https://github.com/shahab-sarmashghi/RESPECT.git
    ```
6. Change to the RESPECT directory and run
    ```
        python setup.py install
    ```

Using RESPECT
-------------
RESPECT accepts uncompressed (*.fastq/.fq/.fa/.fna/.fasta*) or gzipped 
(*.fastq.gz/.fq.gz/.fa.gz/.fna.gz/.fasta.gz*) sequence files 
as well as pre-computed k-mer histograms (`jellyfish histo` *.hist* output
or similarly formatted file) as input. If you have gzipped sequence files
but do not have `gzip` 1.6 or later on your system, use `--decomp` option
to select `zlib` or `gzip` python libraries for decompression (`zlib` is 
preferred method for large inputs).

## Input
### FASTQ files
Files with **.fastq** or **.fq** extension (and **.gz** versions) are 
considered to be a collection of short reads (genome-skims). The k-mer 
histogram is computed for them, and an iterative optimization is run to 
estimate their genomic parameters.

### FASTA files
Files with **.fa**, **.fna**, or **.fasta** extension (and **.gz** versions)
are first processed to determine if they are genome-skims or assemblies 
based on the maximum contig length. Similar to FASTQ files, the k-mer 
histogram is computed for them. If the input is an assembly, the parameters 
are directly computed from the input sequence and the k-mer histogram.

### Histogram files
Files with **.hist** extension can also be provided as the input. The 
file format should be similar to `jellyfish histo` output, a text file 
with two column where the first column is the k-mer multiplicity and the 
second one is the corresponding count of k-mers at that multiplicity:
```
1 20427099
2 216746
3 54799
4 19206
...
```
if you have any histogram file in among you input, you need to provide a
histogram info file using `-I` option. This text file tells RESPECT if the
provided histogram is from a genome-skim or an assembly, and if it is from
a genome-skim, what is the average read length. For example, if you look at
`data\hist_info.txt` you see:
```
Input	read_length
Micromonas_pusilla.hist	-1
Micromonas_pusilla_cov_0.5_err_0.01.hist	100
```
The read length of `Micromonas_pusilla.hist` is `-1`,
a sentinel value to say that the histogram is from an assembly. If a
histogram input information is missing from this file, the input is skipped.

### Input options
For the convenience of user, the input files can be provided in multiple
ways. Using `-i` option, the input files can be provided directly. You can
also provide parent directories using `-d` option and RESPECT takes all 
files with matching extensions as inputs. There is also a `-f` option to
provide a file where each line is a path to an input file. These options 
*are not mutually exclusive* and you can use them simultaneously.

### Mapping file
You can provide a name mapping file through `-m` option to make RESPECT 
use those names in the output. For example, if you open
`data/name_mapping.txt`:
```
Input	Output
Micromonas_pusilla_cov_0.5_err_0.01.fq.gz	mp_fq
GCF_000151265.2_Micromonas_pusilla_CCMP1545_v2.0_genomic.fna	mp_fa
Micromonas_pusilla.hist	mp_ha
Micromonas_pusilla_cov_0.5_err_0.01.hist	mp_hq
```
The first column specifies the input file names (*extensions included*)
and the second column contains the alternative names that should be used
to name the samples in the results. If an input file is missing from this
file, the input name is used for the outputs as well.


## Ouput
The final estimates of parameters are written to `estimated_parameters.txt`
and `estimated_spectra.txt` in the current working directory (you can use 
`-o` option to modify the output directory). The estimated parameters file
headings are:
```
sample	input_type	sequence_type	coverage	genome_length	uniqueness_ratio	HCRM	sequencing_error_rate   average_read_length
``` 
The `input_type` is either `sequence` (for FASTQ/A files) or `histogram` 
(for .hist files). The sequence type determines whether the input is from
an `assembly` or a `genome-skim`. 

The `uniqueness_ratio` is the estimated
ratio of number of unique k-mers (k-mers with one copy) in the *genome* 
to the total length of genome, and is a measure of repetitiveness of the
genome sequence. Low `uniqueness_ratio` values (<0.8) can be a signature of recent
whole genome duplication. 

`HCRM` stands for high copy repeats per million, and is the average
count (per million base-pairs) of the 10 most highly
repetitive k-mers. For assemblies we compute that directly from the k-mer
histogram, and for genome-skims it is computed using the estimated coverage.
High `HCRM` values can be attributed to the presence of transposable elements. 
Please note that the estimate of `HCRM` from genome-skim can be up to 2 times 
the estimate from the assembly because we do not know which strand (forward/reverse) 
of DNA the reads are coming from, so the counts of a k-mer and its reverse 
complement are aggregated.

The `average_read_length` for assemblies is set to `-1` to be consistent with the
convention used in histogram information file (`data\hist_info.txt`).

The estimated repeat spectra are written to `estimated_spectra.txt`. The
column `rX` specifies the estimated number of k-mers with `X` copy in the
genome. You can modify the number of spectra shown in the output using 
`-s` option.

The intermediate results, optimization logs, and k-mer histograms of inputs
are by default written to current_working_directory/tmp, but you can use 
`--tmp` option to specify a non-default path.

## Algorithm options
We use an iterative optimization algorithm to find a set of parameters for
the genome that can explain the observed histogram of k-mers. The size of
k-mers is `31` by default, but you can use smaller k-mers with `-k` option.
Please note that the genome looks more repetitive when using smaller
k-mers, so at some point you might also want to lower the minimum 
permissible uniqueness ratio (`-r` option).

The number of iterations of the algorithm is `1000` by default and can be 
changed using `-N` option. However, be aware that changing that can impact
the speed and convergence of algorithm, so it should be done advisedly.

We do not use a parametric model for the repeat spectrum, and the number 
of spectra can be set by `-n` option, which is `50` by default. Modifying 
this option also impacts the speed of the algorithm as it determines the 
size of the constrained linear program solved in each iteration of the 
algorithm.

## Usage example
After successfully installing RESPECT, change directory to `RESPECT` and
run:
```
respect -d data/ -m data/name_mapping.txt -I data/hist_info.txt -N 10 --debug
```
This command tells RESPECT to process all matching files under `data/`
directory. `data/name_mapping.txt` is used for name mapping, and the
required information about histogram files are provided in 
`data/hist_info.txt`. The algorithm is run for `10` iterations (just for
the purpose of this example) and `--debug` is used to instruct RESPECT to
print out stack traces of any errors happened. The output files will be
written to the current working directory. Run `RESPECT --help` to see all 
available options.
 
[1]: https://www.biorxiv.org/content/10.1101/2021.01.28.428636v1
[2]: http://www.genome.umd.edu/jellyfish.html
[3]: https://www.gurobi.com
[4]: https://conda.io/miniconda.html
[5]: https://github.com/lh3/seqtk
[6]: https://www.gurobi.com/documentation/9.1/quickstart_mac/cs_anaconda_and_grb_conda_.html
[7]: https://www.gurobi.com/documentation/9.1/quickstart_mac/obtaining_a_grb_license.html
