"""Functions and classes to profile k-mer histogram using Jellyfish"""


import os
from subprocess import call, check_output, STDOUT
from shutil import copyfile

from respect import settings


class ProfilerError(Exception):
    """A custom exception used to report errors in k-mer profiling"""


def parse_jellyfish_histo_output(histo_stderr):
    """Parses the output from running jellyfish histo output.

    Parameters
    ----------
    histo_stderr : str
        The jellyfish histo output or similarly formatted histogram

    Returns
    -------
    tuple : (dict, int)
        - count_dict (dict): A dictionary that contains non-zero values of the histogram
        - kmer_sum (int): The total number of k-mers
    """

    count_dict = {0: 0}
    kmer_sum = 0
    for item in histo_stderr.split('\n'):
        if not item.strip():
            continue
        multiplicity = int(float(item.split()[0]))
        count = int(float(item.split()[1]))
        count_dict[multiplicity] = count
        kmer_sum += count * multiplicity

    return count_dict, kmer_sum


def kmer_profiler(input_file, sequence_type, out_name, out_dir, kmer_length, n_threads):
    """Runs Jellyfish on input and returns k-mer histogram and sum.

    Jellyfish count is run on sequence file, then Jellyfish histo is
    called to get the k-mer histogram file. Then,
    parse_jellyfish_histo_output() is called to parse Jellyfish histo
    results and return the histogram dictionary and total sum of kmers.

    Parameters
    ----------
    input_file : str
        The input sequence file
    sequence_type : str
        The type of sequence; 'assembly' or 'genome-skim'
    out_name : str
        The name used for the output files
    out_dir : str
        The directory to write output files
    kmer_length : int
        The size of k-mers
    n_threads : int
        The number of threads used when running jellyfish count

    Raises
    ------
    ProfilerError
        If sequence type is not set yet.
    """

    histo_file = os.path.join(out_dir, out_name + '.hist')
    highest_histo_multiplicity = settings.HIGHEST_JELLYFISH_HISTO_MULTIPLICITY
    mercnt = os.path.join(out_dir, out_name + '.jf')

    # Whether or not counting canonical k-mers based on the sequence type
    if sequence_type == 'assembly':
        call(["jellyfish", "count", "-m", str(kmer_length), "-s", "100M", "-t", str(n_threads), "-o", mercnt,
              input_file], stderr=open(os.devnull, 'w'))
    elif sequence_type == 'genome-skim':
        call(["jellyfish", "count", "-m", str(kmer_length), "-s", "100M", "-t", str(n_threads), "-C", "-o", mercnt,
              input_file], stderr=open(os.devnull, 'w'))
    else:
        raise ProfilerError("The sequence type is not set properly")

    histo_stderr = check_output(["jellyfish", "histo", "-h", str(highest_histo_multiplicity), mercnt], stderr=STDOUT,
                                universal_newlines=True)
    with open(histo_file, mode='w') as f:
        f.write(histo_stderr)
    os.remove(mercnt)

    return parse_jellyfish_histo_output(histo_stderr)


def profile_reader(input_file, out_name, out_dir):
    """Reads input histogram file and returns k-mer histogram and sum"""

    with open(input_file) as f:
        histo_stderr = f.read()

    histo_file = os.path.join(out_dir, out_name + '.hist')
    copyfile(input_file, histo_file)

    return parse_jellyfish_histo_output(histo_stderr)
