import os
import re
import importlib.resources
from importlib import import_module
import configparser
from datetime import timedelta
from dateutil.parser import parse as dateparse
from osgeo import gdal
from cesard.config import (keyval_check, validate_options, validate_value,
                           version_dict as cesard_version_dict)


def get_keys(section):
    """
    get all allowed configuration keys for a section

    Parameters
    ----------
    section: str
        the configuration section to get the allowed keys for.
        Either 'processing' or the name of a SAR processor plugin e.g. 'snap'.

    Returns
    -------
    list[str]
        a list of keys
    """
    if section == 'processing':
        return ['acq_mode', 'annotation', 'aoi_geometry', 'aoi_tiles',
                'ard_dir', 'date_strict', 'db_file', 'dem_type',
                'gdal_threads', 'logfile', 'maxdate',
                'mindate', 'mode', 'processor',
                'sar_dir', 'scene', 'scene_dir', 'sensor',
                'tmp_dir', 'wbm_dir', 'work_dir']
    else:
        try:
            module = import_module(f'asard.{section}')
        except ModuleNotFoundError:
            raise RuntimeError(f"unknown section: {section}.")
        try:
            return module.get_config_keys()
        except AttributeError:
            raise RuntimeError(f"missing function asard.{section}.get_config_keys().")


def read_config_file(config_file=None):
    """
    Reads a configuration file and returns a ConfigParser object

    Parameters
    ----------
    config_file: str or None
        the configuration file name. If None, the default configuration file
        within the package will be used.

    Returns
    -------
    configparser.ConfigParser
    """
    parser = configparser.ConfigParser(allow_no_value=True,
                                       converters={'_datetime': _parse_datetime,
                                                   '_list': _parse_list})
    
    if config_file:
        if not os.path.isfile(config_file):
            raise FileNotFoundError(f"Config file {config_file} does not exist.")
    else:
        with importlib.resources.path(package='asard.resources',
                                      resource='config.ini') as path:
            config_file = str(path)
    
    parser.read(config_file)
    return parser


def get_config(config_file=None, **kwargs):
    """
    Returns the content of a `config.ini` file as a dictionary.

    Parameters
    ----------
    config_file: str or None
        Full path to the config file that should be parsed to a dictionary.
    kwargs: dict[str]
        further keyword arguments overriding configuration found in the config file.

    Returns
    -------
    dict
        Dictionary of the parsed config parameters.
        The keys correspond to the config sections in lowercase letters.
    """
    parser = read_config_file(config_file)
    
    kwargs_proc = {k: v for k, v in kwargs.items() if k in get_keys('processing')}
    
    out = {'processing': _get_config_processing(parser, **kwargs_proc)}
    
    processor_name = out['processing']['processor']
    processor = import_module(f'asard.{processor_name}')
    kwargs_sar = {k: v for k, v in kwargs.items() if k in get_keys(processor_name)}
    out[processor_name] = processor.get_config_section(parser, **kwargs_sar)
    
    return out


def _get_config_processing(parser, **kwargs):
    allowed_keys = get_keys(section='processing')
    try:
        proc_sec = parser['PROCESSING']
    except KeyError:
        msg = "Section 'PROCESSING' does not exist in the config file"
        raise KeyError(msg)
    
    # override config file parameters with additional keyword arguments
    for k, v in kwargs.items():
        if k in allowed_keys:
            proc_sec[k] = v.strip()
    
    # make all relevant paths absolute
    for k in ['work_dir', 'scene', 'scene_dir']:
        v = proc_sec[k]
        proc_sec[k] = 'None' if v in ['', 'None'] else os.path.abspath(v)
    
    # set some defaults
    processing_defaults = {
        'sar_dir': 'SAR',
        'tmp_dir': 'TMP',
        'ard_dir': 'ARD',
        'wbm_dir': 'WBM',
        'gdal_threads': '4',
        'dem_type': 'Copernicus 30m Global DEM',
        'date_strict': 'True',
        'annotation': 'dm,ei,id,lc,li,np,ratio',
        'logfile': 'None'
    }
    
    processing_options = {
        'acq_mode': ['APP', 'APS', 'IMM', 'IMP', 'IMS', 'WSM', 'WSS'],
        'annotation': ['dm', 'ei', 'em', 'id', 'lc',
                       'ld', 'li', 'np', 'ratio', 'wm'],
        'dem_type': ['Copernicus 10m EEA DEM',
                     'Copernicus 30m Global DEM',
                     'Copernicus 30m Global DEM II',
                     'GETASSE30'],
        'mode': ['sar', 'nrb'],
        'sensor': ['ERS1', 'ERS2', 'ASAR']}
    
    for k, v in processing_defaults.items():
        if k not in proc_sec.keys():
            proc_sec[k] = v
    
    # check completeness of configuration parameters
    missing = []
    exclude = ['aoi_tiles', 'aoi_geometry']
    for key in get_keys(section='processing'):
        if key not in proc_sec.keys() and key not in exclude:
            missing.append(key)
    if len(missing) > 0:
        missing_str = '\n - ' + '\n - '.join(missing)
        raise RuntimeError(f"missing the following parameters:{missing_str}")
    
    out = {}
    for k, v in proc_sec.items():
        # check if key is allowed and convert 'None|none|' strings to None
        v = keyval_check(key=k, val=v, allowed_keys=allowed_keys)
        
        if k in ['annotation', 'aoi_tiles', 'mode']:
            v = proc_sec.get_list(k)
        
        validate_value(k, v)
        
        if k == 'mindate' and v is not None:
            v = proc_sec.get_datetime(k)
        if k == 'maxdate' and v is not None:
            date_short = re.search('^[0-9-]{10}$', v) is not None
            v = proc_sec.get_datetime(k)
            if date_short:
                v += timedelta(days=1, microseconds=-1)
        dir_ignore = ['work_dir']
        if k == 'scene_dir' and v is None:
            dir_ignore.append(k)
        if k.endswith('_dir') and k not in dir_ignore:
            if os.path.isabs(v):
                msg = f"Parameter '{k}': '{v}' must be an existing directory"
                assert v is not None and os.path.isdir(v), msg
            else:
                v = os.path.join(proc_sec['work_dir'], v)
        if k.endswith('_file') and not k.startswith('db'):
            msg = f"Parameter '{k}': file {v} could not be found"
            if os.path.isabs(v):
                assert os.path.isfile(v), msg
            else:
                v = os.path.join(proc_sec['work_dir'], v)
                assert os.path.isfile(v), msg
        if k in ['db_file', 'logfile'] and v is not None:
            if not os.path.isabs(v):
                v = os.path.join(proc_sec['work_dir'], v)
        if k == 'gdal_threads':
            v = int(v)
        if k in ['date_strict', 'date_strict']:
            v = proc_sec.getboolean(k)
        
        validate_options(k, v, options=processing_options)
        out[k] = v
    
    # check that a valid scene search option is set
    db_file_set = out['db_file'] is not None
    
    options_set = sum([db_file_set])
    
    if options_set == 0:
        raise RuntimeError("Please define a scene search option.")
    elif options_set > 1:
        raise RuntimeError("Multiple scene search options have been defined. Please choose only one.")
    
    return out


def _parse_datetime(s):
    """Custom converter for configparser:
    https://docs.python.org/3/library/configparser.html#customizing-parser-behaviour"""
    return dateparse(s)


def _parse_list(s):
    """Custom converter for configparser:
    https://docs.python.org/3/library/configparser.html#customizing-parser-behaviour"""
    if s in ['', 'None']:
        return None
    else:
        return [x.strip() for x in s.split(',')]


def gdal_conf(config):
    """
    Stores GDAL configuration options for the current process.
    
    Parameters
    ----------
    config: dict
        Dictionary of the parsed config parameters for the current process.
    
    Returns
    -------
    dict
        Dictionary containing GDAL configuration options for the current process.
    """
    threads = config['processing']['gdal_threads']
    threads_before = gdal.GetConfigOption('GDAL_NUM_THREADS')
    if not isinstance(threads, int):
        raise TypeError("'threads' must be of type int")
    if threads == 1:
        multithread = False
    elif threads > 1:
        multithread = True
        gdal.SetConfigOption('GDAL_NUM_THREADS', str(threads))
    else:
        raise ValueError("'threads' must be >= 1")
    
    return {'threads': threads, 'threads_before': threads_before,
            'multithread': multithread}


def version_dict(processor_name: str) -> dict[str, str]:
    """
    Get the versions of used packages

    Parameters
    ----------
    processor_name:
        the name of the used SAR processor (e.g. `snap`)
    
    Returns
    -------
        a dictionary containing the versions of relevant python packages and
        the return of the `version_dict` function of the used SAR processor.
    """
    import asard
    out = cesard_version_dict()
    out['asard'] = asard.__version__
    processor = import_module(f'asard.{processor_name}')
    out.update(processor.version_dict())
    return out
