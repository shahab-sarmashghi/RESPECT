"""Contains utility functions used by other classes"""


import os
from subprocess import call, check_output, STDOUT
import numpy


def compute_sequence_stat(input_file, out_name, tmp_dir, sample_size=10000):
    """Sub-samples the input sequence file and computes sequence statistics.

    Parameters
    ----------
    input_file : str
        The input FASTQ/A file
    out_name : str
        The name of sub-sampled output file
    tmp_dir : str
        The temporary directory used to write sub-sampled file
    sample_size : int
        The size of sub-sampled sequence file. No sampling if sample_size=-1

    Returns
    -------
    tuple : (float, int, int, int)
        - average_read_length (float): The average length of reads/contigs
        - max_contig_length (int): The longest contig/read length
        - total_length (int): The total length of input sequences
        - total_length_ungapped (int): The total length of input sequences excluding non-ACTG bases
    """

    # out_name = os.path.basename(input_file).rsplit('.f', 1)[0]
    # file_format = os.path.splitext(input_file)[1]
    file_format = '.f' + os.path.basename(input_file).rsplit('.f', 1)[1]

    if sample_size == -1:
        comp_stdout = check_output(["seqtk", "comp", input_file], stderr=STDOUT, universal_newlines=True)
        reads_stat = comp_stdout.split('\n')
    else:
        subsampled_sequence_file = os.path.join(tmp_dir, out_name + '_subsampled' + os.path.splitext(file_format)[0])
        seed = 243
        cmd = "seqtk sample -s {0} {1} {2} > {3}".format(seed, input_file, sample_size, subsampled_sequence_file)
        call(cmd, shell=True)
        comp_stdout = check_output(["seqtk", "comp", subsampled_sequence_file], stderr=STDOUT, universal_newlines=True)
        reads_stat = comp_stdout.split('\n')
        os.remove(subsampled_sequence_file)

    total_length = 0
    total_length_ungapped = 0
    number_of_reads = 0
    max_contig_length = 0
    for stat in reads_stat:
        if not stat.strip():
            continue
        read_length = int(stat.split('\t')[1])
        read_length_ungapped = sum([int(x) for x in stat.split('\t')[2:6]])
        total_length += read_length
        total_length_ungapped += read_length_ungapped
        max_contig_length = max(max_contig_length, read_length)
        number_of_reads += 1

    average_read_length = 1.0 * total_length_ungapped / number_of_reads

    return average_read_length, max_contig_length, total_length, total_length_ungapped


def get_zero_padded_histogram(histogram, max_multiplicity):
    """Extracts and truncates the input histogram.

    Pads the missing entries of input histogram with zero and
    truncates it to max_multiplicity.

    Parameters
    ----------
    histogram : dict
        The compact dictionary of input histogram
    max_multiplicity : int
        The histogram truncation cutoff
    """

    observed_multiplicity_set = set(histogram)
    zero_padded_histogram = []
    for i in range(1, max_multiplicity + 1):
        if i in observed_multiplicity_set:
            zero_padded_histogram.append(float(histogram[i]))
        else:
            zero_padded_histogram.append(0.0)

    return numpy.array(zero_padded_histogram)


def get_uniqueness_ratio(repeat_spectra, genome_length):
    """Returns the uniqueness ratio of repeat spectra of genome.

    If the number of unique k-mers is larger than the genome length,
    the output would be maxed out at 1.0

    Parameters
    ----------
    repeat_spectra : numpy.ndarray
        The estimated repeat spectra
    genome_length : int
        The length of genome
    """

    return min(1.0, 1.0 * repeat_spectra[0] / genome_length)
