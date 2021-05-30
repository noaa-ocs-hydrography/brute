import configparser
import pprint
import io
import glob
import pathlib
import logging
import os
import shutil
import sys
from typing import Union

"""  Sample usage -------

LOGGER = get_logger('xipe.csar.proc')
CONFIG_SECTION = 'combined_raster_processing'

# default_config_name = "default.config"

if __name__ == '__main__':
    if len(sys.argv) > 1:
        use_configs = sys.argv[1:]
    else:
        use_configs = pathlib.Path(__file__).parent.resolve()  # (os.path.dirname(os.path.abspath(__file__))

    warnings = ""
    # LOGGER.info(f'running {len(config_filenames)} configuration(s)')
    for config_filename, config_file in iter_configs(use_configs):
        stringio_warnings = set_stream_logging("xipe", file_level=logging.WARNING, remove_other_file_loggers=False)
        LOGGER.info(f'***************************** Start Run  *****************************')
        LOGGER.info(f'reading "{config_filename}"')
        log_config(config_file, LOGGER)

        config = config_file[CONFIG_SECTION if CONFIG_SECTION in config_file else 'DEFAULT']
"""

def iter_configs(config_filenames:Union[list, str, os.PathLike], log_files:bool=True, default_config_name:Union[str, os.PathLike]=""):
    """ Read all the configs using configparser and optionally modified by a default config file and base_configs.
    A ConfigParser object is created.  Then the default_config_name is loaded, if applicable.
    Then loads all configs from the base_configs directory (local to the script) listed in the [DEFAULT] section 'additional_configs' entry.
    Then looks for a subdirectory of the logged in user (os.getlogin()) and uses that if it exists, otherwise uses the current dir.
    Finally iterates each config in the user subdirectory or current directory so they have the highest priority.

    Parameters
    ----------
    config_filenames
        Either an iterable list of config paths or a str/PathLike that points to a directory to scan.
        If only one config is desired then send it as a list ["c:\\test\\example.config"]
    log_files
        If True then a pair of log files will be created in a log directory in the same path as the config.
        One log gets all message levels the other is Warnings and Errors.
    default_config_name
        The default name to load with each config.  Any info in the current config will take priority over the default values.

    Returns
    -------
    filename str, ConfigParser instance

    """
    user_directory_exists = False
    if isinstance(config_filenames, (str, os.PathLike)):
        use_configs = pathlib.Path(config_filenames)  # (os.path.dirname(os.path.abspath(__file__))
        user_dir = use_configs.joinpath(os.getlogin())
        if user_dir.exists():
            user_directory_exists = True
            use_configs = user_dir
        config_filenames = use_configs.glob('*.config')
    # remove the default names from the list of files to process
    config_filenames = [pathlib.Path(p) for p in config_filenames]
    for config_filename in filter(lambda fname: fname.name != default_config_name, config_filenames):
        config_path, just_cfg_filename = config_filename.parent, config_filename.name
        if user_directory_exists:
            base_config_path = config_path.parent.joinpath("base_configs")  # parallel to the user directory, so up one directory from config
        else:
            base_config_path = config_path.joinpath("base_configs")  # local to the config file
        os.makedirs(os.path.join(config_path, "logs"), exist_ok=True)
        if log_files:
            set_file_logging("nbs", os.path.join(config_path, 'logs',
                                                  just_cfg_filename + ".log"))  # sets the parent logger to output to a file
            set_file_logging("nbs", os.path.join(config_path, 'logs', just_cfg_filename + ".warnings.log"), file_level=logging.WARNING,
                             remove_other_file_loggers=False)  # sets the parent logger to output to a file
        config_file = configparser.ConfigParser(interpolation=configparser.ExtendedInterpolation())
        config_file.read(config_filename)
        if default_config_name:
            config_file.read(os.path.join(config_path, default_config_name))
        try:
            extra_confs = [fname.strip() for fname in config_file['DEFAULT']['additional_configs'].split(",")]
            extra_confs.reverse()  # file should be in most significant to least - read in the opposite order
        except KeyError:
            extra_confs = []
        for fname in extra_confs:
            if os.sep not in fname and "/" not in fname and "\\" not in fname:
                fname = base_config_path.joinpath(fname)
            config_file.read(fname)
        config_file.read(config_filename)  # read the original file again so it's defaults are used and not overwritten by the sub-configs

        yield(config_filename, config_file)

def log_config(config_file, LOGGER, absolute=True):
    """ Writes the data from a config file into a log - so the parameters used are stored in the processing log.

    Parameters
    ----------
    config_file
        configparser.ConfigParser instance
    LOGGER
        logging.Logger instance
    absolute
        flag for if variables are evaluated and written with values or as variables like ${var_name}

    Returns
    -------
    None

    """
    if absolute:
        ss = pprint.pformat({section: dict(config_file[section]) for section in config_file.sections()}, width=500)
        LOGGER.warning(ss)
    else:  # show with variables in the data e.g.  ${section:variable}
        # archive the config used
        ss = io.StringIO()
        config_file.write(ss)
        ss.seek(0)
        LOGGER.info(ss.read())
        del ss



def get_logger(name: str, log_filename: str = None, file_level: int = None, console_level: int = None, log_format: str = None) -> logging.Logger:
    if console_level is None:
        console_level = logging.INFO
    logger = logging.getLogger(name)

    # check if logger is already configured
    if logger.level == logging.NOTSET and len(logger.handlers) == 0:
        # check if logger has a parent
        if '.' in name:
            logger.parent = get_logger(name.rsplit('.', 1)[0])
        else:
            # otherwise create a new split-console logger
            logger.setLevel(logging.DEBUG)
            if console_level != logging.NOTSET:
                # this creates a logger which writes DEBUG and INFO to stdout
                # and another that writes WARNING, ERROR, CRITICAL to stderr
                # stderr gets colored differently in pycharm and could be read from a pipe if needed
                # otherwise it doesn't have much effect
                if console_level <= logging.INFO:
                    class LoggingOutputFilter(logging.Filter):
                        def filter(self, rec):
                            return rec.levelno in (logging.DEBUG, logging.INFO)

                    console_output = logging.StreamHandler(sys.stdout)
                    console_output.setLevel(console_level)
                    console_output.addFilter(LoggingOutputFilter())
                    logger.addHandler(console_output)

                console_errors = logging.StreamHandler(sys.stderr)
                console_errors.setLevel(max((console_level, logging.WARNING)))
                logger.addHandler(console_errors)
            # Will create a log file on disk, if a filename is specified
            set_file_logging(name, log_filename, file_level)

    if log_format is None:
        log_format = '[%(asctime)s] %(name)-15s %(levelname)-8s: %(message)s'
    log_formatter = logging.Formatter(log_format)
    for handler in logger.handlers:
        handler.setFormatter(log_formatter)

    return logger


def set_stream_logging(logger_name: str, file_level: int = None, log_format: str = None,
                     remove_other_file_loggers: bool = True):
    io_string = io.StringIO()
    log_stream = logging.StreamHandler(io_string)
    set_file_logging(logger_name, log_stream, file_level, log_format, remove_other_file_loggers)
    return io_string

def set_file_logging(logger_name: str, log_file: Union[str, logging.StreamHandler] = None, file_level: int = None, log_format: str = None,
                     remove_other_file_loggers: bool = True):
    logger = logging.getLogger(logger_name)
    logger.info(f"{logger.name} saving to {log_file}")
    if remove_other_file_loggers:  # remove string and file loggers except for the default stderr, stdout
        for existing_file_handler in [handler for handler in logger.handlers if isinstance(handler, (logging.FileHandler, logging.StreamHandler))]:
            if existing_file_handler.stream not in (sys.stderr, sys.stdout):
                logger.removeHandler(existing_file_handler)
    if log_file is not None:
        if isinstance(log_file, str):
            file_handler = logging.FileHandler(log_file)
        else:
            file_handler = log_file
        if file_level is None:
            file_level = logging.DEBUG
        file_handler.setLevel(file_level)
        if log_format is None:
            log_format = '[%(asctime)s] %(name)-15s %(levelname)-8s: %(message)s'
        log_formatter = logging.Formatter(log_format)
        file_handler.setFormatter(log_formatter)
        logger.addHandler(file_handler)
        pass
