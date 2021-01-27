#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import os
import multiprocessing
import sys

from respect import settings
from respect import __version__
from respect.respect_functions import run_respect


def main():
    # Arguments parser
    parser = argparse.ArgumentParser(description='Estimating the genomic parameters for input genome-skims, '
                                                 'or pre-computed histograms. Input files can be provided '
                                                 'directly, by their parent directory, and file of files (see '
                                                 'input options). Results are by default written to the current '
                                                 'working directory into estimated-parameters.txt and '
                                                 'estimated-spectra.txt as tab-separated tables. Use -o option to '
                                                 'change the output directory', add_help=False)

    # Input/output options
    parser_io_group = parser.add_argument_group("Input/output options", description='Different input options can be '
                                                                                    'used simultaneously')
    parser_io_group.add_argument('-i', '--input', nargs='+', type=str, metavar='<files>',
                                 help='Input sequences (FASTQ/FASTA) or histogram (.hist) files. To '
                                      'provide multiple files, separate them with a space: respect -i input1 '
                                      'input2 ...')

    parser_io_group.add_argument('-d', '--input-directories', nargs='+', type=str, metavar='<paths>',
                                 help='Directory of input sequences (FASTQ/FASTA) or histogram (.hist) files. Every '
                                      'file under the input directory (not recursive) is examined and a warning is '
                                      'issued if the extension is not valid. To provide multiple directories, '
                                      'separate them with a space: respect -d dir1 dir2 ...')

    parser_io_group.add_argument('-f', '--input-file', type=str, metavar='<file>',
                                 help='List of files. A file where each line is the path to an input file')

    parser_io_group.add_argument('-m', '--mapping-file', type=str, metavar='<file>',
                                 help='Tab-separated table where the first column is the name of input files ('
                                      'including file extensions like a.fq or b.hist) and the second column specifies '
                                      'the output name. The first line/row should be column names/headers. For the '
                                      'samples not included in the mapping, input names are used')

    parser_io_group.add_argument('-I', '--info-file', type=str, metavar='<file>',
                                 help='Tab-separated table where the first column is the names of histogram (.hist) '
                                      'files and the second column specifies the average read length of the original '
                                      'genome-skim. If the histogram is computed from an assembly, the second column '
                                      'should be -1 so respect knows the sequence type. The first line/row should be '
                                      'column names/headers. The input .hist files not included in this file will be '
                                      'ignored')

    parser_io_group.add_argument('-o', '--output-directory', default=os.getcwd(), type=str, metavar='<path>',
                                 help='Directory of output files. Default: current working directory')

    parser_io_group.add_argument('-s', '--spectra-output-size', default=settings.SPECTRA_OUTPUT_SIZE,
                                 type=int, metavar='<int>',
                                 help='The number of spectra written to the output file. Default: {}. The last bin '
                                      'collects every k-mer with higher multiplicity as well (Please note that this '
                                      'option is only for the output and different from the number of spectra used in '
                                      'the algorithm which is set by -n option)'.format(settings.SPECTRA_OUTPUT_SIZE))

    # Algorithm options
    parser_algorithm_group = parser.add_argument_group("Algorithm options", description='Please note that changing '
                                                                                        'the default parameters '
                                                                                        'can impact the accuracy and '
                                                                                        'speed of the algorithm')
    parser_algorithm_group.add_argument('-k', '--kmer-size', default=settings.KMER_SIZE, type=int,
                                        choices=list(range(1, 32)), metavar='<int>',
                                        help='K-mer length. It should be an integer between 1 and 31. Default: {}. '
                                             'You might need to lower the minimum permissible uniqueness ratio (see -r '
                                             'option) if you use small k-mer sizes'.format(settings.KMER_SIZE))

    parser_algorithm_group.add_argument('-N', '--iterations', default=settings.ITERATIONS, type=int, metavar='<int>',
                                        help='Number of iterations of the optimization algorithm. '
                                             'Default: {}'.format(settings.ITERATIONS))

    parser_algorithm_group.add_argument('-T', '--temperature', default=settings.TEMPERATURE_COEFFICIENT, type=float,
                                        metavar='<float>', help='Temperature coefficient of the simulated annealing '
                                                                'scheme used for the optimization. '
                                                                'Default: {}'.format(settings.TEMPERATURE_COEFFICIENT))

    parser_algorithm_group.add_argument('-r', '--min-r1l', default=settings.MINIMUM_R1L, type=float,
                                        metavar='<float>', help='Smallest permissible uniqueness ratio (r_1/L). '
                                                                'Default: {}'.format(settings.MINIMUM_R1L))

    parser_algorithm_group.add_argument('-l', '--norm', default=settings.ERROR_FUNCTION_NORM, type=float,
                                        metavar='<float>', help='Norm used in computing the error function between '
                                                                'observed and expected k-mer counts. '
                                                                'Default: {}'.format(settings.ERROR_FUNCTION_NORM))

    parser_algorithm_group.add_argument('-n', '--spectra-number', default=settings.SPECTRA_NUMBER, type=int,
                                        metavar='<int>', help='Highest multiplicity considered in modeling the repeat '
                                                              'spectrum. Default: {}'.format(settings.SPECTRA_NUMBER))

    # Other options
    parser_other_group = parser.add_argument_group("Other options")
    parser_other_group.add_argument('--threads', default=multiprocessing.cpu_count(), type=int,
                                    metavar='<int>', choices=list(range(1, multiprocessing.cpu_count() + 1)),
                                    help='Maximum number of processing cores to use. '
                                         'Default for this machine: {0}'.format(multiprocessing.cpu_count()))

    # parser_other_group.add_argument('--log', default=os.path.join(os.getcwd(), 'respect.log'), type=str,
    #                                 metavar='<file>', help='Log file. Default: current_working_directory/respect.log')

    parser_other_group.add_argument('--log', type=str, metavar='<file>', help='Log file. Default: printed to the '
                                                                              'standard error if no file is provided')

    parser_other_group.add_argument('--tmp', default=os.path.join(os.getcwd(), 'tmp'), type=str, metavar='<path>',
                                    help='Parent directory to store optimization logs, k-mer histograms, and estimated '
                                         'parameters during the convergence of algorithm. A separate directory is '
                                         'created per each sample under this parent directory. Default: '
                                         'current_working_directory/tmp')

    parser_other_group.add_argument('--version', action='version', version='{} {}'.format(parser.prog, __version__),
                                    help='Print the version number')

    parser_other_group.add_argument('--debug', action='store_true',
                                    help='Print the traceback when an exception is raised')

    parser_other_group.add_argument('-h', '--help', action='help', help='show this help message and exit')

    args = parser.parse_args()

    # Handling traceback on exceptions
    def exception_handler(exception_type, exception, traceback, debug_hook=sys.excepthook):
        if args.debug:
            debug_hook(exception_type, exception, traceback)
        else:
            print("{0}: {1}".format(exception_type.__name__, exception))

    sys.excepthook = exception_handler

    run_respect(args)


if __name__ == '__main__':
    main()
