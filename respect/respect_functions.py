import os
import errno
import fnmatch
import logging
import pandas
import multiprocessing

from respect import settings
from respect.paramter_estimator import ParameterEstimator


def get_parameter_estimator_results(parameter_estimator, n_spectra, parameters_dataframe, spectra_dataframe):
    """Extracts the estimated parameters from estimator instance.

    Returns data-frames of estimated parameters and repeat spectra with
    the results of input ParameterEstimator instance appended to them.

    Parameters
    ----------
    parameter_estimator : ParameterEstimator
        The ParameterEstimator instance to get the estimates from
    n_spectra : int
        The number of spectra shown in the output
    parameters_dataframe : pandas.DataFrame
        The dataframe to store estimated parameters of input files
    spectra_dataframe : pandas.DataFrame
        The dataframe to store estimated repeat spectra of input files
    """

    coverage = parameter_estimator.coverage
    r1l = parameter_estimator.uniqueness_ratio
    hcrm = parameter_estimator.high_copy_repeats_per_million
    e = parameter_estimator.sequencing_error_rate

    parameters_dataframe = parameters_dataframe.append(
        {'sample': parameter_estimator.output_name, 'input_type': parameter_estimator.input_type,
         'sequence_type': parameter_estimator.sequence_type,
         'coverage': format(coverage, "0.2f") if coverage else coverage,
         'genome_length': parameter_estimator.genome_length,
         'uniqueness_ratio': format(r1l, ".2f") if r1l else r1l,
         'HCRM': format(hcrm, ".2f") if hcrm else hcrm,
         'sequencing_error_rate': format(e, ".4f") if e else e}, ignore_index=True)

    # Making the last bin collect all k-mers with n_spectra and higher multiplicity
    compact_repeat_spectra = [int(parameter_estimator.repeat_spectra[i]) for i in range(n_spectra - 1)]
    compact_repeat_spectra += [int(sum(parameter_estimator.repeat_spectra[n_spectra - 1:]))]

    spectra_dataframe = spectra_dataframe.append(pandas.Series([parameter_estimator.output_name] +
                                                               compact_repeat_spectra, index=spectra_dataframe.columns),
                                                 ignore_index=True)

    return parameters_dataframe, spectra_dataframe


def get_histograms_info(input_files, histogram_info_file):
    """Reads and processes the histogram info file.

    Returns a dictionary of read length information and a list of missing
    .hist files. For duplicate entries in the file, warning is issued and
    the first instance is used.

    Parameters
    ----------
    input_files : list
        The list of input files, each is a string
    histogram_info_file : str
        The file which provides read length and sequence type information of
        input histogram files
    """

    hist_info_dict = {}
    missing_hist_files = []

    if histogram_info_file:
        hist_info_dataframe = pandas.read_csv(histogram_info_file, sep='\t')
        info_series = pandas.Series(hist_info_dataframe.iloc[:, 1].values, index=hist_info_dataframe.iloc[:, 0])

        if info_series.index.unique().__len__() < info_series.index.__len__():
            logging.warning("Duplicate input names in the histogram info file; the first instance is used")
        info_series_dedup = info_series[~info_series.index.duplicated()]

        for in_file in input_files:
            if in_file.endswith('.hist'):
                in_name = os.path.basename(in_file)
                in_name_no_ext = os.path.splitext(in_name)[0]
                if in_name in info_series_dedup.index:
                    hist_info_dict[in_file] = info_series_dedup[in_name]
                elif in_name_no_ext in info_series_dedup.index:
                    hist_info_dict[in_file] = info_series_dedup[in_name_no_ext]
                else:
                    logging.warning("{0} is missing from the histogram info file; it's skipped".format(in_file))
                    missing_hist_files.append(in_file)

    elif any([in_file.endswith('.hist') for in_file in input_files]):
        logging.warning("The histogram info file is not provided, please see -I option. All .hist files are skipped")
        missing_hist_files = [in_file for in_file in input_files if in_file.endswith('.hist')]

    return hist_info_dict, missing_hist_files


def get_output_files_names(input_files, mapping_file):
    """Returns a dictionary which maps input files to output names.

    For duplicate entries in the mapping, warning is issued and the first
    is used.

    Parameters
    ----------
    input_files : list
        The list of input files, each is a string
    mapping_file : str
        The file which specifies the output names for all or some of inputs
    """

    # Initializing the deduplicated mapping series and output name dictionary
    mapping_final = pandas.Series()
    output_names = {}

    if mapping_file:
        name_mapping_dataframe = pandas.read_csv(mapping_file, sep='\t')
        mapping_series = pandas.Series(name_mapping_dataframe.iloc[:, 1].values,
                                       index=name_mapping_dataframe.iloc[:, 0])

        if mapping_series.index.unique().__len__() < mapping_series.index.__len__():
            logging.warning('Duplicate input names in the mapping file; the first instance is used')
        mapping_series_index_dedup = mapping_series[~mapping_series.index.duplicated()]

        if mapping_series_index_dedup.unique().__len__() < mapping_series_index_dedup.__len__():
            logging.warning('Duplicate output names in the mapping file; the first instance is used')
        mapping_final = mapping_series_index_dedup.drop_duplicates()

    for in_file in input_files:
        in_name = os.path.basename(in_file)
        if in_name in mapping_final.index:
            output_names[in_file] = mapping_final[in_name]
        else:
            output_names[in_file] = in_name

    return output_names


def custom_formatwarning(msg, *args, **kwargs):
    """Monkeypatch the format of warning message.

    Parameters
    ----------
    msg : str
        The custom message
    """

    # ignore everything except the message
    return 'Warning: ' + str(msg) + '\n'


def valid_file(in_file):
    """Returns True if the input is a file and has valid extension.

    Warning is issued otherwise. Please note that it only checks file
    extension and does not validate the content or format of input.

    Parameters
    ----------
    in_file : str
        The input file
    """

    if not os.path.isfile(in_file):
        logging.warning("{} is not a file; it's skipped".format(in_file))
        return False
    elif True not in (fnmatch.fnmatch(in_file.rsplit('.gz', 1)[0], '*' + ext) for ext in
                      settings.VALID_FILE_EXTENSIONS):
        logging.warning("{} does not have valid extension; it's skipped".format(in_file))
        return False
    return True


def run_respect(args):
    """Processes the input arguments and sets up RESPECT algorithm.

    Parameters
    ----------
    args : namespace
        The namespace object returned by the parse_args() method

    Raises
    ------
    SystemExit
        If no valid input file is provided.
    """

    if all(x is None for x in [args.input, args.input_directories, args.input_file]):
        raise SystemExit("No input specified. Please provide an input FASTQ/A or .hist file")

    # Setting up logging
    log_format = '%(asctime)s %(levelname)s:%(message)s'
    if args.log:
        print("Using Gurobi Optimizer to solve linear programs. License information:")
        logging.basicConfig(format=log_format, filename=args.log, level=logging.DEBUG, filemode='w')
    else:
        logging.basicConfig(format=log_format, level=logging.DEBUG, filemode='w')

    # Making a list of input files provided through different input options
    input_sequence_files = []

    # -i option; files directly provided to the input
    if args.input:
        input_sequence_files += [in_file for in_file in args.input if valid_file(in_file)]

    # -d option; files under the provided directories
    if args.input_directories:
        input_sequence_files += [os.path.join(in_dir, f) for in_dir in args.input_directories
                                 for f in os.listdir(in_dir) if valid_file(os.path.join(in_dir, f))]

    # -f option; files from the list of files
    if args.input_file:
        with open(args.input_file) as file_files:
            input_sequence_files += [line.strip() for line in file_files if valid_file(line.strip())]

    if not input_sequence_files:
        raise SystemExit("No valid input file is provided. Please provide an input FASTQ/A or .hist file")

    unique_input_files = list(set([os.path.abspath(f) for f in input_sequence_files]))

    # -m option; mapping to output names
    output_names = get_output_files_names(unique_input_files, args.mapping_file)

    # -I option; histogram info file
    histograms_info, missing_info_histograms = get_histograms_info(unique_input_files, args.info_file)
    filtered_input_files = list(set(unique_input_files).difference(missing_info_histograms))

    # --tmp option; making a parent directory for temporary files
    try:
        os.makedirs(args.tmp)
    except OSError as Error:
        if Error.errno != errno.EEXIST:
            raise

    # Initiating estimator instances for each input file
    parameter_estimators_list = []
    for input_file in filtered_input_files:
        logging.info("Processing {}...".format(output_names[input_file]))
        try:
            parameter_estimator = ParameterEstimator(input_file, args.kmer_size, output_names[input_file], args.tmp)
            parameter_estimator.preprocess_input_file(histograms_info)
            parameter_estimator.set_kmer_histogram(args.threads)
        except:
            logging.exception("Error occurred when processing {}; it's skipped".format(input_file))
        else:
            parameter_estimators_list.append(parameter_estimator)

    # Initiating output data-frames
    parameters_dataframe = pandas.DataFrame(columns=['sample', 'input_type', 'sequence_type', 'coverage',
                                                     'genome_length', 'uniqueness_ratio', 'HCRM',
                                                     'sequencing_error_rate'])
    spectra_dataframe = pandas.DataFrame(columns=['sample'] +
                                                 ['r{}'.format(j) for j in range(1, args.spectra_output_size + 1)])

    # Estimating parameters for preprocessed files in parallel
    number_of_worker_processes = min(args.threads, parameter_estimators_list.__len__())
    pool = multiprocessing.Pool(number_of_worker_processes)
    results = [pool.apply_async(parameter_estimator,
                                args=(args.spectra_number, args.norm, args.iterations,
                                      args.min_r1l, args.temperature))
               for parameter_estimator in parameter_estimators_list]

    for result in results:
        try:
            parameter_estimator = result.get(9999999)
        except:
            logging.error("Error occurred while trying to get estimated parameters for a sample")
        else:
            parameters_dataframe, spectra_dataframe = get_parameter_estimator_results(parameter_estimator,
                                                                                      args.spectra_output_size,
                                                                                      parameters_dataframe,
                                                                                      spectra_dataframe)
    pool.close()
    pool.join()

    # Writing the results to the output files
    logging.info("Writing the results to the output files...")
    parameters_dataframe.to_csv(os.path.join(args.output_directory, 'estimated-parameters.txt'), sep='\t', mode='w',
                                index=False, na_rep="NA")
    spectra_dataframe.to_csv(os.path.join(args.output_directory, 'estimated-spectra.txt'), sep='\t', mode='w',
                             index=False, na_rep="NA")
