"""Defines classes used to structure the parameter estimation"""


import os
import errno
import logging

from respect.timer import Timer
from respect.profiling import kmer_profiler, profile_reader
from respect.utils import compute_sequence_stat
from respect import settings
from respect.optimizer import Optimizer
from respect.utils import get_zero_padded_histogram, get_uniqueness_ratio


class EstimatorError(Exception):
    """A custom exception used to report errors in estimating parameters"""


class ParameterEstimator:
    """
    A class used to represent input sequence/histogram files

    ...

    Attributes
    ----------
    _histogram : dict
        The histogram of observed counts
    _total_kmers : int
        The total number of observed k-mers
    input_type : str
        The type of input file: "histogram" or "sequence"
    sequence_type : str
        The sequence type of input file: "assembly" or "genome_skim"
    read_length : float
        The average length of sequencing reads
    genome_length : int
        The estimated length of genome
    genome_length_ungapped : int
        The length of genome excluding gaps (N's) when input type is assembly
    coverage : float
        The estimated coverage of genome-skim
    sequencing_error_rate : float
        The estimated sequencing error rate of genome-skim
    repeat_spectra : numpy.ndarray
        The estimated repeat spectra of input
    uniqueness_ratio : float
        The ratio of unique k-mers to genome length
    high_copy_repeats_per_million : float
        The average count (per million base-pairs) of the 10 most highly
        repetitive k-mers
    input_file : str
        The input file
    tmp_dir : str
        The directory to save temporary files and output logs, as well
        as the iteration results
    kmer_length : int
        The size of k-mers
    output_name : str
        The name to use for the corresponding output files

    Methods
    -------
    preprocess_input_file(histograms_info)
        Processes the input file and sets its basic attributes
    compute_kmer_histogram(n_threads)
        Computes the k-mer histogram for sequence inputs
    read_kmer_histogram()
        Reads the k-mer histogram for histogram inputs
    set_kmer_histogram(n_threads)
        Sets the histogram by calling respective functions
    compute_assembly_parameters()
        Computes the assembly parameters
    estimate_genome_skim_parameters(spectra_number, error_norm,
                                    iterations_number, min_r1l,
                                    temperature)
        Estimates the genome-skim parameters
    estimate_genomic_parameters(spectra_number, error_norm,
                                iterations_number, min_r1l,
                                temperature)
        Estimates the parameters by calling respective functions
    """

    def __init__(self, input_file, kmer_length, output_name, parent_tmp_dir):
        self._histogram = None
        self._total_kmers = None
        self.input_type = None
        self.sequence_type = None
        self.read_length = None
        self.genome_length = None
        self.genome_length_ungapped = None
        self.coverage = None
        self.sequencing_error_rate = None
        self.repeat_spectra = None
        self.uniqueness_ratio = None
        self.high_copy_repeats_per_million = None
        self.input_file = input_file
        self.kmer_length = kmer_length
        self.output_name = output_name
        self.tmp_dir = os.path.join(parent_tmp_dir, output_name)

        # making a directory for temporary files
        try:
            os.makedirs(self.tmp_dir)
        except OSError as Error:
            if Error.errno != errno.EEXIST:
                raise

    # @Timer(logger=logging.info)
    def preprocess_input_file(self, histograms_info):
        """Processes the input file and sets its basic attributes.

        Sets the input type based on the extension of input file. For
        input histograms, the sequence type and read length (if any)
        are set from the input histograms_info dictionary. For sequence
        files, the sequence type is determined based on the maximum
        length of any contiguous sequence inside the input. The read
        length and genome length are set depending on the sequence
        type.

        Parameters
        ----------
        histograms_info : dict
            The dictionary of read length values for .hist files
        """

        if self.input_file.endswith('.hist'):
            self.input_type = 'histogram'
            # need to get hist content checked
            read_length = histograms_info[self.input_file]
            if read_length == -1:
                self.sequence_type = 'assembly'
            else:
                self.sequence_type = 'genome-skim'
                self.read_length = read_length

        else:
            self.input_type = 'sequence'
            if self.input_file.endswith('a'):
                stats = compute_sequence_stat(self.input_file, self.output_name, self.tmp_dir, sample_size=-1)
                average_read_length = stats[0]
                maximum_contig_length = stats[1]
                total_sequence_length = stats[2]
                total_sequence_length_ungapped = stats[3]
                if maximum_contig_length > settings.ASSEMBLY_CONTIG_LENGTH_THRESHOLD:
                    self.sequence_type = 'assembly'
                    self.genome_length = total_sequence_length
                    self.genome_length_ungapped = total_sequence_length_ungapped
                else:
                    self.sequence_type = 'genome-skim'
                    self.read_length = average_read_length
            else:
                stats = compute_sequence_stat(self.input_file, self.output_name, self.tmp_dir)
                average_read_length = stats[0]
                self.sequence_type = 'genome-skim'
                self.read_length = average_read_length

    @Timer(logger=logging.info)
    def compute_kmer_histogram(self, n_threads):
        """Computes the k-mer histogram for sequence inputs.

        Sets _histogram and _total_kmers attributes.

        Parameters
        ----------
        n_threads : int
            The number of threads used when profiling sequence inputs
        """

        profiler_output = kmer_profiler(self.input_file, self.sequence_type, self.output_name, self.tmp_dir,
                                        self.kmer_length, n_threads)
        self._histogram = profiler_output[0]
        self._total_kmers = profiler_output[1]

    # @Timer(logger=logging.info)
    def read_kmer_histogram(self):
        """Reads the k-mer histogram for histogram inputs.

        Sets _histogram and _total_kmers attributes. In addition
        if the histogram is from an assembly, the genome length is
        set too.
        """

        reader_output = profile_reader(self.input_file, self.output_name, self.tmp_dir)
        self._histogram = reader_output[0]
        self._total_kmers = reader_output[1]
        if self.sequence_type == 'assembly':
            self.genome_length = reader_output[1]
            self.genome_length_ungapped = reader_output[1]

    def set_kmer_histogram(self, n_threads):
        """Sets the histogram by calling respective functions.

        Depending on the input type, the histogram is computed (for
        sequence input) or read from input (for hist input).

        Parameters
        ----------
        n_threads : int
            The number of threads passed to compute_kmer_histogram()

        Raises
        ------
        EstimatorError
            If input type is not set yet.
        """

        if self.input_type == 'histogram':
            self.read_kmer_histogram()
        elif self.input_type == 'sequence':
            self.compute_kmer_histogram(n_threads)
        else:
            raise EstimatorError("The input type is not set properly")

    # @Timer(logger=logging.info)
    def compute_assembly_parameters(self, spectra_number):
        """Computes the assembly parameters.

        Sets the repeat spectra directly from the _histogram. Also,
        uniqueness ratio and HCRM values are computed accordingly.

        Parameters
        ----------
        spectra_number : int
            The number of spectra considered for the repeat spectra

        Raises
        ------
        EstimatorError
            If the assembly length is not set yet.
        """

        if self.genome_length is None:
            raise EstimatorError("The assembly length is not set yet")

        self.repeat_spectra = get_zero_padded_histogram(self._histogram, spectra_number)
        self.uniqueness_ratio = get_uniqueness_ratio(self.repeat_spectra, self.genome_length)
        multiplicity_list = sorted(list(self._histogram))
        high_copy_repeats_average_count = 1.0 * sum(multiplicity_list[-settings.HCRM_NUMBER:]) / settings.HCRM_NUMBER
        self.high_copy_repeats_per_million = 10**6 * high_copy_repeats_average_count / self.genome_length

    @Timer(logger=logging.info)
    def estimate_genome_skim_parameters(self, spectra_number, error_norm, iterations_number, min_r1l, temperature):
        """Estimates the genome-skim parameters.

        D

        Parameters
        ----------
        spectra_number : int
            The number of spectra considered for the repeat spectra
        error_norm : float
            The norm of error function used in the iterative estimation
        iterations_number : int
            The number of iterations of the iterative algorithm
        min_r1l : float
            The minimum permissible uniqueness ratio of the estimated
            repeat spectra
        temperature : float
            The temperature coefficient of the iterative algorithm

        Raises
        ------
        EstimatorError
            If read length is not set yet.
        """

        if self.read_length is None:
            raise EstimatorError("The read length is not set yet")

        optimizer = Optimizer(self._histogram, self.kmer_length, self._total_kmers, self.read_length, spectra_number,
                              error_norm, self.output_name, self.tmp_dir)
        optimizer.set_initial_values()
        optimizer.run_simulated_annealing(iterations_number, min_r1l, temperature)

        self.genome_length = optimizer.genome_length
        self.coverage = optimizer.coverage
        self.sequencing_error_rate = optimizer.sequencing_error_rate
        self.repeat_spectra = optimizer.repeat_spectra
        self.uniqueness_ratio = optimizer.uniqueness_ratio

        multiplicity_list = sorted(list(self._histogram))
        high_copy_repeats_average_count = 1.0 * sum(multiplicity_list[-settings.HCRM_NUMBER:]) / settings.HCRM_NUMBER
        high_copy_repeats_average_count_corrected = high_copy_repeats_average_count / optimizer.error_free_kmer_coverage
        self.high_copy_repeats_per_million = 10**6 * high_copy_repeats_average_count_corrected / self.genome_length

    def estimate_genomic_parameters(self, spectra_number, error_norm, iterations_number, min_r1l, temperature):
        """Estimates the parameters by calling respective functions.

        Depending on the sequence type, the parameters are either
        computed (for assemblies) or estimated using the iterative
        algorithm (for genome-skims)

        Parameters
        ----------
        spectra_number : int
            The number of spectra considered for the repeat spectra
        error_norm : float
            The norm of error function used in the iterative estimation
        iterations_number : int
            The number of iterations of the estimation algorithm
        min_r1l : float
            The minimum allowed uniqueness ratio of estimated spectra
        temperature : float
            The temperature coefficient of the iterative algorithm

        Raises
        ------
        EstimatorError
            If _histogram or sequence type are not set yet.
        """

        if self._histogram is None:
            raise EstimatorError("The histogram is not set yet")

        if self.sequence_type == 'assembly':
            self.compute_assembly_parameters(spectra_number)
        elif self.sequence_type == 'genome-skim':
            self.estimate_genome_skim_parameters(spectra_number, error_norm, iterations_number, min_r1l, temperature)
        else:
            raise EstimatorError("The sequence type is not set properly")

        return self

    def __call__(self, *args, **kwargs):
        """To make the instance callable and therefore pickleable.

        Multiprocessing needs the function passed to worker processes
        to be pickleable and instance methods are. So, instead of
        passing self.estimate_genomic_parameters(), we pass the
        ParameterEstimator class instance itself. Also, error handling
        added here to provide informative stack traces when using pool
        of workers and an exception is raised in any worker process.
        """

        try:
            return self.estimate_genomic_parameters(*args, **kwargs)
        except Exception as error:
            logging.exception("Error occurred when estimating parameters for {}; it's skipped".format(
                self.input_file))
            raise error
