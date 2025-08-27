import os
import sys
import logging
from textwrap import dedent
from osgeo import gdal
import spatialist
import pyroSAR
from pyroSAR import examine
import ERS_NRB

log = logging.getLogger('s1ard')


def set_logging(config, debug=False):
    """
    Set logging for the current process.

    Parameters
    ----------
    config: dict
        Dictionary of the parsed config parameters for the current process.
    debug: bool
        Set logging level to DEBUG?

    Returns
    -------
    logging.Logger
        The log handler for the current process.
    """
    level = logging.DEBUG if debug else logging.INFO
    
    logger = logging.getLogger('ERS_NRB')
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
    
    # add s1ard logger
    log_s1ard = logging.getLogger('s1ard')
    log_s1ard.setLevel(level)
    log_s1ard.addHandler(handler)
    
    return logger


def _log_process_config(logger, config):
    """
    Adds a header to the logfile, which includes information about the current processing configuration.

    Parameters
    ----------
    logger: logging.Logger
        The logger to which the header is added to.
    config: dict
        Dictionary of the parsed config parameters for the current process.
    """
    try:
        snap_config = examine.ExamineSnap()
        core = snap_config.get_version('core')
        microwavetbx = snap_config.get_version('microwavetbx')
        snap_core = f"{core['version']} | {core['date']}"
        snap_microwavetbx = f"{microwavetbx['version']} | {microwavetbx['date']}"
    except RuntimeError:
        snap_core = 'unknown'
        snap_microwavetbx = 'unknown'
    
    header = f"""
    ====================================================================================================================
    PROCESSING CONFIGURATION

    mode                {config['processing']['mode']}
    aoi_tiles           {config['processing']['aoi_tiles']}
    aoi_geometry        {config['processing']['aoi_geometry']}
    mindate             {config['processing']['mindate'].isoformat()}
    maxdate             {config['processing']['maxdate'].isoformat()}
    acq_mode            {config['processing']['acq_mode']}

    work_dir            {config['processing']['work_dir']}
    sar_dir             {config['processing']['sar_dir']}
    tmp_dir             {config['processing']['tmp_dir']}
    wbm_dir             {config['processing']['wbm_dir']}
    
    scene_dir           {config['processing']['scene_dir']}
    db_file             {config['processing']['db_file']}
    dem_type            {config['processing']['dem_type']}
    gdal_threads        {config['processing']['gdal_threads']}

    ====================================================================================================================
    SOFTWARE

    ERS_NRB             {ERS_NRB.__version__}
    snap-core           {snap_core}
    snap-microwavetbx   {snap_microwavetbx}
    python              {sys.version}
    python-pyroSAR      {pyroSAR.__version__}
    python-spatialist   {spatialist.__version__}
    python-GDAL         {gdal.__version__}

    ====================================================================================================================
    """
    logger.info(dedent(header))
