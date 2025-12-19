import os
import sys
import logging
from typing import Any
from datetime import datetime
from asard.config import version_dict

log = logging.getLogger('asard')


def set_logging(
        config: dict[str, Any],
        debug: bool = False
) -> logging.Logger:
    """
    Set logging for the current process.

    Parameters
    ----------
    config:
        Dictionary of the parsed config parameters for the current process.
    debug:
        Set logging level to DEBUG?

    Returns
    -------
        The log handler for the current process.
    """
    level = logging.DEBUG if debug else logging.INFO
    
    logger = logging.getLogger('asard')
    logger.setLevel(level)
    
    log_format = "[%(asctime)s] [%(levelname)5s] %(message)s"
    formatter = logging.Formatter(fmt=log_format,
                                  datefmt='%Y-%m-%d %H:%M:%S')
    
    logfile = config['processing']['logfile']
    if logfile is not None:
        os.makedirs(os.path.dirname(logfile), exist_ok=True)
        handler = logging.FileHandler(filename=logfile, mode='a')
    else:
        handler = logging.StreamHandler(sys.stdout)
    logger.addHandler(handler)
    
    # Add header first with simple formatting
    formatter_simple = logging.Formatter("%(message)s")
    handler.setFormatter(formatter_simple)
    _log_process_config(logger=logger, config=config)
    
    # Use normal formatting from here on
    handler.setFormatter(formatter)
    
    # add pyroSAR logger
    log_pyro = logging.getLogger('pyroSAR')
    log_pyro.setLevel(level)
    log_pyro.addHandler(handler)
    
    # add cesard logger
    log_cesard = logging.getLogger('cesard')
    log_cesard.setLevel(level)
    log_cesard.addHandler(handler)
    
    return logger


def _log_process_config(
        logger: logging.Logger,
        config: dict[str, Any]
) -> None:
    """
    Adds a header to the logfile, which includes information about
    the current processing configuration.

    Parameters
    ----------
    logger:
        The logger to which the header is added to.
    config:
        Dictionary of the parsed config parameters for the current process.
    """
    processor_name = config['processing']['processor']
    
    sw_versions = version_dict(processor_name=processor_name)
    
    max_len_sw = len(max(sw_versions.keys(), key=len))
    max_len_main = len(max(config['processing'].keys(), key=len))
    max_len_proc = len(max(config[processor_name].keys(), key=len))
    max_len = max(max_len_sw, max_len_main, max_len_proc) + 4
    
    lines = []
    lines.append('=' * 100)
    for section in ['PROCESSING', processor_name.upper()]:
        lines.append(f'{section}')
        for k, v in config[section.lower()].items():
            if k == 'dem_prepare_mode':
                continue
            if isinstance(v, datetime):
                val = v.isoformat()
            else:
                val = v
            lines.append(f"{k: <{max_len}}{val}")
        lines.append('=' * 100)
    lines.append('SOFTWARE')
    for k, v in sw_versions.items():
        lines.append(f"{k: <{max_len}}{v}")
    lines.append('=' * 100)
    header = '\n'.join(lines)
    logger.info(header)
