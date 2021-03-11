import logging
import sys
import os.path as path

def_log_file_level=logging.INFO
def_log_stream_level=logging.INFO
def_log_file_path=None

formatter = logging.Formatter('%(asctime)s - %(name)s - %(module)s - %(funcName)s - %(levelname)s - %(message)s')

def initialize_logger(logger_identifier, log_file_path=def_log_file_path, log_stream_level=def_log_stream_level,
                      log_file_level=def_log_file_level):
    logger = logging.getLogger(logger_identifier)
    logger.setLevel(logging.DEBUG)
    if log_file_path is None:
        filepath = f'{logger_identifier}.log'
    else:
        filepath = path.join(log_file_path, f'{logger_identifier}.log')

    fh = logging.FileHandler(filepath)
    fh.setLevel(log_file_level)
    # create console handler with a higher log level
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(log_stream_level)
    # create formatter and add it to the handlers
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)
    # add the handlers to the logger
    logger.addHandler(fh)
    logger.addHandler(ch)

    return logger


def getlogger(logger_identifier):
    return logging.getLogger(logger_identifier)

def set_stream_level(level):
    ch.setLevel(level)
    
def set_logfile_level(level):
    fh.setLevel(log_file_level)

def set_logfile_path(path):
    log_file_path = path
    if log_file_path is None:
        filepath = logger_identifier
    else:
        filepath = path.join(log_file_path, f'{logger_identifier}.log')
    fh = logging.FileHandler('filepath')
    
    
