"""Defines classes used to structure optimization"""

import numpy as np
from scipy.stats import poisson
from numpy.linalg import norm
import pandas as pd
import os
import gurobipy
import logging

from respect import settings
from respect.utils import get_zero_padded_histogram, get_uniqueness_ratio


class OptimizerError(Exception):
    """A custom exception used to report errors in estimating parameters"""


class Optimizer:
    """
    A class used to represent the optimization algorithm

    ...

    Attributes
    ----------
    _truncated_histogram : numpy.ndarray
        The histogram truncated based on the initial coverage
    _multiplicity_cutoff : int
        The histogram truncation cutoff; number of observed counts
    _error_function : float
        The objective function of optimization
    _parameters_dataframe : pandas.DataFrame
        The dataframe to store the parameters during iterations
    _spectra_dataframe : pandas.DataFrame
        The dataframe to store the repeat spectra during iterations
    error_free_kmer_coverage : float
        The estimated error free k-mer coverage
    genome_length : int
        The estimated length of genome
    coverage : float
        The estimated base coverage
    sequencing_error_rate : float
        The estimated sequencing error rate
    repeat_spectra : numpy.ndarray
        The estimated repeat spectra
    uniqueness_ratio : float
        The ratio of unique k-mers to genome length
    histogram : dict
        The histogram of observed k-mer counts
    kmer_length : int
        The size of k-mers
    total_kmers : int
        The total number of k-mers
    read_length : float
        The average length of input reads
    spectra_number : int
        The number of spectra considered in the repeat spectra
    error_norm : float
        The norm used in computing the error function
    output_name : str
        The name to use for the corresponding output files
    tmp_dir : str
        The directory to save output files

    Methods
    -------
    get_poisson_matrix(error_free_kmer_coverage)
        Returns the poisson coefficient matrix
    get_initial_kmer_coverage(index)
        Returns the initial estimate of k-mer coverage
    get_error_function(repeat_spectra, poisson_matrix)
        Returns the error function of optimization
    set_initial_values()
        Setting the initial values of parameters
    apply_error_correction()
        Applies the sequencing error rate correction
    save_iteration_results(iteration)
        Saves the results of each iteration into data-frames
    solve_constrained_lp_using_gurobi(observed_counts, coefficients_matrix, low, up)
        Solves for the repeat spectra using Gurobi LP solver
    estimate_repeat_spectra(histogram, poisson_matrix)
        Returns the estimated repeat spectra
    run_simulated_annealing(iterations_number, min_r1l, temperature)
        Runs the simulated annealing algorithm to minimize the error function
    """

    def __init__(self, histogram, kmer_length, total_kmers, read_length, spectra_number, error_norm, output_name,
                 tmp_dir):
        self._truncated_histogram = None
        self._multiplicity_cutoff = None
        self._error_function = None
        self._parameters_dataframe = None
        self._spectra_dataframe = None
        self.error_free_kmer_coverage = None
        self.genome_length = None
        self.coverage = None
        self.sequencing_error_rate = None
        self.repeat_spectra = None
        self.uniqueness_ratio = None
        self.histogram = histogram
        self.kmer_length = kmer_length
        self.total_kmers = total_kmers
        self.read_length = read_length
        self.spectra_number = spectra_number
        self.error_norm = error_norm
        self.output_name = output_name
        self.tmp_dir = tmp_dir

    def get_poisson_matrix(self, error_free_kmer_coverage):
        """Returns the poisson coefficient matrix.

        Parameters
        ----------
        error_free_kmer_coverage : float
            The error-free k-mer coverage
        """

        poisson_matrix = np.zeros((self._multiplicity_cutoff, self.spectra_number))
        for j in range(1, self.spectra_number + 1):
            for i in range(1, self._multiplicity_cutoff + 1):
                poisson_matrix[i - 1][j - 1] = poisson.pmf(i, j * error_free_kmer_coverage)
        return poisson_matrix

    def get_initial_kmer_coverage(self, index):
        """Returns the initial estimate of k-mer coverage.

        Parameters
        ----------
        index : int
            The maximizing index of the histogram
        """
        x = self.error_free_kmer_coverage
        o = self._truncated_histogram
        return (np.exp(-x) * (x ** (index + 1)) / np.math.factorial(index + 1)) * o[0] / o[index] + x * (1 - np.exp(-x))

    def get_error_function(self, repeat_spectra, poisson_matrix):
        """Returns the error function of optimization.

        The error function is set to the norm of difference between
        observed and modeled k-mer counts.

        Parameters
        ----------
        repeat_spectra : numpy.ndarray
            The estimated repeat spectra
        poisson_matrix : numpy.ndarray
            The poisson coefficients matrix
        """

        p_times_r = poisson_matrix.dot(np.array(repeat_spectra))
        return norm(self._truncated_histogram[1:] - p_times_r[1:], ord=self.error_norm)

    def set_initial_values(self):
        """Setting the initial values of parameters.

        The initial estimate of coverage and other parameters are
        computed based on the assumption that genome is not repetitive.
        """

        initial_padded_histogram_counts = get_zero_padded_histogram(self.histogram, settings.MAX_COVERAGE)
        most_frequent_multiplicity_index = 1 + initial_padded_histogram_counts[1:].argmax()
        index = most_frequent_multiplicity_index

        try:
            self.error_free_kmer_coverage = 1.0 * (index + 2) * initial_padded_histogram_counts[index + 1] / \
                                            initial_padded_histogram_counts[index]
        except IndexError:
            raise OptimizerError("The histogram index is out of range; cannot estimate the initial coverage")

        self._multiplicity_cutoff = int(self.spectra_number * max(1, self.error_free_kmer_coverage))
        self._truncated_histogram = get_zero_padded_histogram(self.histogram, self._multiplicity_cutoff)

        poisson_matrix = self.get_poisson_matrix(self.error_free_kmer_coverage)
        kmer_coverage = self.get_initial_kmer_coverage(index)

        self.genome_length = int(self.total_kmers / kmer_coverage)

        erroneous_kmers_count = self.genome_length * (kmer_coverage - self.error_free_kmer_coverage)

        self.repeat_spectra = np.zeros(self.spectra_number)
        self.repeat_spectra[0] = self.genome_length
        self.uniqueness_ratio = 1.0

        self.coverage = kmer_coverage * self.read_length / (self.read_length - self.kmer_length + 1)

        error_free_kmer_prob = self.error_free_kmer_coverage / kmer_coverage
        self.sequencing_error_rate = max(0, 1 - error_free_kmer_prob ** (1.0 / self.kmer_length))

        self._error_function = self.get_error_function(self.repeat_spectra, poisson_matrix)

    def apply_error_correction(self):
        """Applies the sequencing error rate correction

        If the estimated coverage is less than the threshold, the
        correction learned from training data is read from file and
        applied to the estimated error rate. The estimates of coverage
        and genome length are updated accordingly.
        """

        if self.coverage < settings.ERROR_CORRECTION_COVERAGE_THRESHOLD:
            training_data_dir = settings.TRAINING_DATA_DIR
            error_correction = pd.read_csv(os.path.join(training_data_dir, 'sequencing_error_rate_correction.txt'),
                                           sep='\t', header=0, index_col=0, float_precision='round_trip')

            ratio = round(np.clip(self.uniqueness_ratio, 1e-4, 1.0), 4)
            corrected_sequencing_error_rate = self.sequencing_error_rate / (1 + error_correction.loc[ratio,
                                                                                                     'relative_error'])

            self.sequencing_error_rate = corrected_sequencing_error_rate
            kmer_coverage = self.error_free_kmer_coverage / ((1 - self.sequencing_error_rate) ** self.kmer_length)
            self.coverage = kmer_coverage * self.read_length / (self.read_length - self.kmer_length + 1)
            self.genome_length = int(self.total_kmers / kmer_coverage)

    def save_iteration_results(self, iteration):
        """Saves the results of each iteration into data-frames.

        The data-frames are also written to files located in tmp_dir
        so the convergence of the algorithm can be inspected if needed.

        Parameters
        ----------
        iteration : int
            The iteration number of the simulated annealing algorithm
        """

        self._parameters_dataframe = self._parameters_dataframe.append(
            {'iteration': iteration,
             'coverage': format(self.coverage, ".2f"),
             'genome_length': int(self.genome_length),
             'uniqueness_ratio': format(self.uniqueness_ratio, ".2f"),
             'error_function': format(self._error_function, ".2f")}, ignore_index=True)

        self._spectra_dataframe = self._spectra_dataframe.append(pd.Series([iteration] + list(self.repeat_spectra),
                                                                           index=self._spectra_dataframe.columns),
                                                                 ignore_index=True)

        self._parameters_dataframe.to_csv(os.path.join(self.tmp_dir, 'estimated_parameters_iterations.txt'),
                                          sep='\t', mode='w', index=False)
        self._spectra_dataframe.to_csv(os.path.join(self.tmp_dir, 'estimated_spectra_iterations.txt'),
                                       sep='\t', mode='w', index=False)

    @staticmethod
    def solve_constrained_lp_using_gurobi(observed_counts, coefficients_matrix, lower, upper, log_file):
        """Solves for the repeat spectra using Gurobi LP solver.

        First, the equations are set based on o=Pr, and lower and
        upper bounds on the ratio of consecutive spectra are set to
        form a linear program. Then, the optimizer is called to solve
        for the repeat spectra.

        Parameters
        ----------
        observed_counts : numpy.ndarray
            The array of observed counts; constants of linear equations
        coefficients_matrix : numpy.ndarray
            The poisson matrix; coefficients of linear equations
        lower : pandas.Series
            The lower bound on the ratio of consecutive spectra
        upper : pandas.Series
            The upper bound on the ratio of consecutive spectra
        log_file : str
            The upper bound on the ratio of consecutive spectra
        """

        o = observed_counts
        P = coefficients_matrix
        try:
            model = gurobipy.Model('Spectra estimation')

            # Sending logs to log_file and suppressing output to console
            # m.setParam('OutputFlag', 0)
            model.setParam('LogToConsole', 0)
            model.setParam('LogFile', log_file)

            (I, H) = P.shape
            r = model.addVars(H, name='r')
            u = model.addVars(I, obj=1.0, name='u')
            model.addConstrs((-u[i] - r.prod({k: v for k, v in enumerate(P[i])}) <= -o[i]
                              for i in range(I)), 'binR')
            model.addConstrs((-u[i] + r.prod({k: v for k, v in enumerate(P[i])}) <= o[i]
                              for i in range(I)), 'binL')

            model.addConstrs((r[j] - upper[j] - upper[j] * r[j + 1] <= 0 for j in range(H - 1)), 'specR')
            model.addConstrs((r[j] - lower[j] - lower[j] * r[j + 1] >= 0 for j in range(H - 1)), 'specL')

            model.optimize()

            if model.status == gurobipy.GRB.Status.OPTIMAL:
                solution = model.getAttr('x', r)
                return np.array(list(solution.values()))
        except gurobipy.GurobiError as e:
            print(e.message)

    def estimate_repeat_spectra(self, histogram, poisson_matrix):
        """Returns the estimated repeat spectra.

        Given the observed histogram and poisson matrix, sets up the
        constraints and runs a linear programming solver to get the
        repeat spectra. Then, the spectral ratios learned using splines
        are used to refine the estimated spectra.

        Parameters
        ----------
        histogram : numpy.ndarray
            The observed k-mer counts
        poisson_matrix : numpy.ndarray
            The poisson coefficients matrix
        """

        training_data_dir = settings.TRAINING_DATA_DIR
        n_splines = settings.SPLINE_NUMBER

        spectral_ratio_lower_bound = pd.Series(pd.read_csv(os.path.join(training_data_dir,
                                                                        'spectral_ratio_lower_bounds.txt'),
                                                           header=None).values.flatten())
        spectral_ratio_upper_bound = pd.Series(pd.read_csv(os.path.join(training_data_dir,
                                                                        'spectral_ratio_upper_bounds.txt'),
                                                           header=None).values.flatten())

        constrained_spectra = self.solve_constrained_lp_using_gurobi(histogram, poisson_matrix,
                                                                     spectral_ratio_lower_bound,
                                                                     spectral_ratio_upper_bound,
                                                                     os.path.join(self.tmp_dir, 'gurobi.log'))

        splines_list = [pd.read_csv(os.path.join(training_data_dir, 'spline_{}.txt'.format(i + 1)), sep='\t', header=0,
                                    index_col=0, float_precision='round_trip') for i in range(n_splines)]

        spectral_residuals = [1.0 * constrained_spectra[i] / norm(constrained_spectra[i:], ord=1) for i in range(
            n_splines)]

        spline_constrained_spectra = constrained_spectra
        for i in range(n_splines):
            residual = round(np.clip(spectral_residuals[i], 1e-4, 1.0), 4)
            spline_constrained_spectra[i + 1] = int(1.0 * spline_constrained_spectra[i] / splines_list[i].loc[residual,
                                                                                                              'ratio'])

        return spline_constrained_spectra

    def run_simulated_annealing(self, iterations_number, min_r1l, temperature):
        """Runs the simulated annealing algorithm to minimize the error function.

        At each iteration, a new candidate coverage is randomly selected within
        a neighbourhood of current estimate, and a repeat spectra is computed
        using Spline Linear Programming method accordingly. Candidate coverages
        are accepted/rejected based on the resulting error function following a
        simulated annealing scheme. The estimated parameters and repeat spectra
        at each iteration are written to files.

        Parameters
        ----------
        iterations_number : int
            The number of iterations of simulated annealing
        min_r1l : float
            The minimum permissible r1/L value
        temperature : float
            The temperature parameter of simulated annealing

        Raises
        ------
        OptimizerError
            If the algorithm is not initialized properly.
        """

        # Initializing the iteration results dataframes
        self._parameters_dataframe = pd.DataFrame(columns=['iteration', 'coverage', 'genome_length', 'uniqueness_ratio',
                                                           'error_function'])
        self._spectra_dataframe = pd.DataFrame(columns=['iteration'] +
                                                       [str(j) for j in range(1, self.spectra_number + 1)])

        # Saving initial estimates/values; iteration 0
        self.save_iteration_results(0)

        cc = settings.CONTRACTION_COEFFICIENT
        random_number_gen_1 = np.random.RandomState(settings.RANDOM_NEIGHBOUR_SEED)
        random_number_gen_2 = np.random.RandomState(settings.ACCEPTANCE_SEED)

        # Defining variables with shorter names for readability
        lambda_ef = self.error_free_kmer_coverage
        epsilon = self.sequencing_error_rate
        o = self._truncated_histogram
        initial_error_function = self._error_function

        logging.info("Starting iterations to estimate parameters of {0}".format(self.output_name))
        for iteration in range(1, iterations_number + 1):
            lambda_ef_next = lambda_ef + cc * random_number_gen_1.uniform(low=-lambda_ef / 2, high=2 * lambda_ef)
            kmer_coverage = lambda_ef_next / ((1 - epsilon) ** self.kmer_length)
            poisson_matrix_next = self.get_poisson_matrix(lambda_ef_next)
            repeat_spectra_next = self.estimate_repeat_spectra(o[1:], poisson_matrix_next[1:, :])

            # Protecting the run against failure in an iteration
            if not np.isnan(repeat_spectra_next).any():
                genome_length_next = int(self.total_kmers / kmer_coverage)
                p_times_r_next = poisson_matrix_next.dot(np.array(repeat_spectra_next))
                erroneous_kmers_count = o[0] - p_times_r_next[0]
                coverage_next = kmer_coverage * self.read_length / (self.read_length - self.kmer_length + 1)

                error_function_next = self.get_error_function(repeat_spectra_next, poisson_matrix_next)
                error_function_next += initial_error_function * np.heaviside(min_r1l - 1.0 * repeat_spectra_next[0]
                                                                             / genome_length_next, 0)
                acceptance_probability = 1
                if error_function_next > self._error_function:
                    acceptance_probability = np.exp(-(error_function_next - self._error_function) / (
                            temperature * iterations_number / iteration))

                if error_function_next < self._error_function or acceptance_probability > random_number_gen_2.uniform():
                    lambda_ef = lambda_ef_next
                    self.error_free_kmer_coverage = lambda_ef_next
                    self._error_function = error_function_next
                    self.coverage = coverage_next
                    self.genome_length = genome_length_next
                    self.repeat_spectra = repeat_spectra_next
                    self.uniqueness_ratio = get_uniqueness_ratio(repeat_spectra_next, genome_length_next)
                    # epsilon = epsilon_next
                else:
                    cc *= settings.CONTRACTION_COEFFICIENT

            # Saving the current state
            self.save_iteration_results(iteration)

        # Applying the correction on the sequencing error
        # The genome length and coverage are modified accordingly
        self.apply_error_correction()
