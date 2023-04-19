#! /usr/bin/env python3


"""Provide functionalities for logging, the command line parser, helper functions 
for interacting with the shell, etc. 
"""


import os
import subprocess
import sys

import logging
logger = logging.getLogger(__name__)


# parser and logging


def add_file_options(parser):

    file_options = parser.add_argument_group('file options')

    file_options.add_argument('-o', '--overwrite', help="""If this flag is present, then 
        existing files will be overwritten.""", action='store_true')

    file_options.add_argument('-k', '--keep', help="""If this flag is present, then 
        intermediary files are not deleted.""", action='store_true')


def add_logging_options(parser, default_log_file=""):
    """ This function add options for logging to an argument parser. In 
        particular, it adds options for logging to a file, stdout and stderr.
        In addition, it adds options for controlling the logging level of each
        of the loggers, and a general option for controlling all of the loggers.

        Args:
            parser (argparse.ArgumentParser): an argument parser

        Returns:
            None, but the parser has the additional options added
    """

    logging_options = parser.add_argument_group("logging options")

    default_log_file = ""
    logging_level_choices = ['NOTSET', 'DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL']
    default_logging_level = 'WARNING'
    default_specific_logging_level = 'NOTSET'

    logging_options.add_argument('--log-file', help="This option specifies a file to "
        "which logging statements will be written (in addition to stdout and "
        "stderr, if specified)", default=default_log_file)
    logging_options.add_argument('--log-stdout', help="If this flag is present, then "
        "logging statements will be written to stdout (in addition to a file "
        "and stderr, if specified)", action='store_true')
    logging_options.add_argument('--no-log-stderr', help="Unless this flag is present, then "
        "logging statements will be written to stderr (in addition to a file "
        "and stdout, if specified)", action='store_true')

    logging_options.add_argument('--logging-level', help="If this value is specified, "
        "then it will be used for all logs", choices=logging_level_choices,
        default=default_logging_level)
    logging_options.add_argument('--file-logging-level', help="The logging level to be "
        "used for the log file, if specified. This option overrides "
        "--logging-level.", choices=logging_level_choices, 
        default=default_specific_logging_level)
    logging_options.add_argument('--stdout-logging-level', help="The logging level to be "
        "used for the stdout log, if specified. This option overrides "
        "--logging-level.", choices=logging_level_choices, 
        default=default_specific_logging_level)
    logging_options.add_argument('--stderr-logging-level', help="The logging level to be "
        "used for the stderr log, if specified. This option overrides "
        "--logging-level.", choices=logging_level_choices, 
        default=default_specific_logging_level)


def update_logging(args, logger=None, 
        format_str='%(levelname)-8s %(name)-8s %(asctime)s : %(message)s'):

    """ This function interprets the logging options in args. Presumably, these
        were added to an argument parser using add_logging_options.

    Parameters
    ----------
    args: argparse.Namespace
        a namespace with the arguments added by add_logging_options

    logger: logging.Logger or None
        a logger which will be updated. If None is given, then the default
        logger will be updated.

    format_str: string
        The logging format string. Please see the python logging documentation
        for examples and more description.

    Returns
    -------
    None, but the default (or given) logger is updated to take into account
        the specified logging options
    """

    # find the root logger if another logger is not specified
    if logger is None:
        logger = logging.getLogger('')
            
    logger.handlers = []

    # set the base logging level
    level = logging.getLevelName(args.logging_level)
    logger.setLevel(level)

    # now, check the specific loggers

    if len(args.log_file) > 0:
        h = logging.FileHandler(args.log_file)
        formatter = logging.Formatter(format_str)
        h.setFormatter(formatter)
        if args.file_logging_level != 'NOTSET':
            l = logging.getLevelName(args.file_logging_level)
            h.setLevel(l)
        logger.addHandler(h)

    if args.log_stdout:
        h = logging.StreamHandler(sys.stdout)
        formatter = logging.Formatter(format_str)
        h.setFormatter(formatter)
        if args.stdout_logging_level != 'NOTSET':
            l = logging.getLevelName(args.stdout_logging_level)
            h.setLevel(l)
        logger.addHandler(h)

    log_stderr = not args.no_log_stderr
    if log_stderr:
        h = logging.StreamHandler(sys.stderr)
        formatter = logging.Formatter(format_str)
        h.setFormatter(formatter)
        if args.stderr_logging_level != 'NOTSET':
            l = logging.getLevelName(args.stderr_logging_level)
            h.setLevel(l)
        logger.addHandler(h)
        
        
# shell


def check_programs_exist(programs, raise_on_error=True, package_name=None, 
            logger=logger):

    """ This function checks that all of the programs in the list cam be
        called from python. After checking all of the programs, an exception
        is raised if any of them are not callable. Optionally, only a warning
        is raised. The name of the package from which the programs are
        available can also be included in the message.

        Internally, this program uses shutil.which, so see the documentation
        for more information about the semantics of calling.

        Arguments:
            programs (list of string): a list of programs to check

        Returns:
            list of string: a list of all programs which are not found

        Raises:
            EnvironmentError: if any programs are not callable, then
                an error is raised listing all uncallable programs.
    """
    import shutil

    missing_programs = []
    for program in programs:
        exe_path = shutil.which(program)

        if exe_path is None:
            missing_programs.append(program)

    if len(missing_programs) > 0:
        missing_programs = ' '.join(missing_programs)
        msg = "The following programs were not found: " + missing_programs

        if package_name is not None:
            msg = msg + ("\nPlease ensure the {} package is installed."
                .format(package_name))

        if raise_on_error:
            raise EnvironmentError(msg)
        else:
            logger.warning(msg)

    return missing_programs


def check_call_step(cmd, current_step = -1, init_step = -1, call=True, 
        raise_on_error=True):
    
    logging.info(cmd)
    ret_code = 0

    if current_step >= init_step:
        if call:
            #logging.info(cmd)
            logging.info("calling")
            ret_code = subprocess.call(cmd, shell=True)

            if raise_on_error and (ret_code != 0):
                raise subprocess.CalledProcessError(ret_code, cmd)
            elif (ret_code != 0):
                msg = ("The command returned a non-zero return code\n\t{}\n\t"
                    "Return code: {}".format(cmd, ret_code))
                logger.warning(msg)
        else:
            msg = "skipping due to --do-not-call flag"
            logging.info(msg)
    else:
        msg = "skipping due to --init-step; {}, {}".format(current_step, init_step)
        logging.info(msg)

    return ret_code



def check_call(cmd, call=True, raise_on_error=True):
    return check_call_step(cmd, call=call, raise_on_error=raise_on_error)


def call_if_not_exists(cmd, out_files, in_files=[], overwrite=False, call=True,
            raise_on_error=True, file_checkers=None, num_attempts=1, 
            to_delete=[], keep_delete_files=False):

    """ This function checks if out_file exists. If it does not, or if overwrite
        is true, then the command is executed, according to the call flag.
        Otherwise, a warning is issued stating that the file already exists
        and that the cmd will be skipped.

        Additionally, a list of input files can be given. If given, they must
        all exist before the call will be executed. Otherwise, a warning is 
        issued and the call is not made.

        However, if call is False, the check for input files is still made,
        but the function will continue rather than quitting. The command will
        be printed to the screen.

        The return code from the called program is returned.

        By default, if the called program returns a non-zero exit code, an
        exception is raised.

        Furthermore, a dictionary can be given which maps from a file name to
        a function which check the integrity of that file. If any of these
        function calls return False, then the relevant file(s) will be deleted
        and the call made again. The number of attempts to succeed is given as
        a parameter to the function.

        Args:
            cmd (string): the command to execute

            out_files (string or list of strings): path to the files whose existence 
                to check. If they do not exist, then the path to them will be 
                created, if necessary.

            in_files (list of strings): paths to files whose existence to check
                before executing the command

            overwrite (bool): whether to overwrite the file (i.e., execute the 
                command, even if the file exists)

            call (bool): whether to call the command, regardless of whether the
                file exists

            raise_on_error (bool): whether to raise an exception on non-zero 
                return codes

            file_checkers (dict-like): a mapping from a file name to a function
                which is used to verify that file. The function should return
                True to indicate the file is okay or False if it is corrupt. The
                functions must also accept "raise_on_error" and "logger" 
                keyword arguments.

            num_attempts (int): the number of times to attempt to create the
                output files such that all of the verifications return True.

            to_delete (list of strings): paths to files to delete if the command
                is executed successfully

            keep_delete_files (bool): if this value is True, then the to_delete
                files will not be deleted, regardless of whether the command
                succeeded

        Returns:
            int: the return code from the called program

        Warnings:
            warnings.warn if the out_file already exists and overwrite is False
            warnings.warn if the in_files do not exist

        Raises:
            subprocess.CalledProcessError: if the called program returns a
                non-zero exit code and raise_on_error is True
                
            OSError: if the maximum number of attempts is exceeded and the 
                file_checkers do not all return true and raise_on_error is True

        Imports:
            os
            shell
    """
    import os
    import shlex

    ret_code = 0

    # check if the input files exist
    missing_in_files = []
    for in_f in in_files:
        # we need to use shlex to ensure that we remove surrounding quotes in
        # case the file name has a space, and we are using the quotes to pass
        # it through shell
        in_f = shlex.split(in_f)[0]

        if not os.path.exists(in_f):
            missing_in_files.append(in_f)

    if len(missing_in_files) > 0:
        msg = "Some input files {} are missing. Skipping call: \n{}".format(missing_in_files, cmd)
        logger.warn(msg)
        return ret_code

        # This is here to create a directory structue using "do_not_call". In
        # hindsight, that does not seem the best way to do this, so it has been
        # removed.
        #if call:
        #    return


    # make sure we are working with a list
    if isinstance(out_files, str):
        out_files = [out_files]

    # check if the output files exist
    all_out_exists = False
    if out_files is not None:
        all_out_exists = all([os.path.exists(of) for of in out_files])

    all_valid = True
    if overwrite or not all_out_exists:
        attempt = 0
        while attempt < num_attempts:
            attempt += 1

            # create necessary paths
            if out_files is not None:
                [os.makedirs(os.path.dirname(x), exist_ok=True) for x in out_files]
            
            # make the call
            ret_code = check_call(cmd, call=call, raise_on_error=raise_on_error)

            # do not check the files if we are not calling anything
            if (not call) or (file_checkers is None):
                break

            # now check the files
            all_valid = True
            for filename, checker_function in file_checkers.items():
                msg = "Checking file for validity: {}".format(filename)
                logger.debug(msg)

                is_valid = checker_function(filename, logger=logger, 
                                raise_on_error=False)

                # if the file is not valid, then rename it
                if not is_valid:
                    all_valid = False
                    invalid_filename = "{}.invalid".format(filename)
                    msg = "Rename invalid file: {} to {}".format(filename, invalid_filename)
                    logger.warning(msg)

                    os.rename(filename, invalid_filename)

            # if they were all valid, then we are done
            if all_valid:
                break


    else:
        msg = "All output files {} already exist. Skipping call: \n{}".format(out_files, cmd)
        logger.warn(msg)

    # now, check if we succeeded in creating the output files
    if not all_valid:
        msg = ("Exceeded maximum number of attempts for cmd. The output files do "
            "not appear to be valid: {}".format(cmd))

        if raise_on_error:
            raise OSError(msg)
        else:
            logger.critical(msg)

    elif (not keep_delete_files):
        # the command succeeded, so delete the specified files
        for filename in to_delete:
            if os.path.exists(filename):
                msg = "Removing file: {}".format(filename)
                logger.info(cmd)
                
                os.remove(filename)

    return ret_code


# utils


def check_keys_exist(d, keys):
    """ This function ensures the given keys are present in the dictionary. It
        does not other validate the type, value, etc., of the keys or their
        values. If a key is not present, a KeyError is raised.

        The motivation behind this function is to verify that a config dictionary
        read in at the beginning of a program contains all of the required values.
        Thus, the program will immediately detect when a required config value is
        not present and quit.

        Input:
            d (dict) : the dictionary

            keys (list) : a list of keys to check
        Returns:
            list of string: a list of all programs which are not found

        Raises:
            KeyError: if any of the keys are not in the dictionary
    """
    missing_keys = [k for k in keys if k not in d]

    
    if len(missing_keys) > 0:
        missing_keys = ' '.join(missing_keys)
        msg = "The following keys were not found: " + missing_keys
        raise KeyError(msg)

    return missing_keys


