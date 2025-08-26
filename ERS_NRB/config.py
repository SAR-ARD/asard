import os
import re
import importlib.resources
from importlib import import_module
import configparser
from datetime import datetime, timedelta
from osgeo import gdal
from s1ard.config import keyval_check, validate_options, validate_value


def get_keys(section):
    """
    get all allowed configuration keys

    Parameters
    ----------
    section: {'processing'}
        the configuration section to get the allowed keys for.

    Returns
    -------
    list[str]
        a list of keys
    """
    if section == 'processing':
        return ['acq_mode', 'annotation', 'aoi_geometry', 'aoi_tiles', 'ard_dir',
                'db_file', 'dem_type', 'gdal_threads', 'logfile',
                'maxdate', 'mindate', 'mode', 'processor', 'sar_dir', 'scene_dir',
                'tmp_dir', 'wbm_dir', 'work_dir']
    else:
        raise RuntimeError(f"unknown section: {section}. Options: 'processing'.")


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
                                       converters={'_annotation': _parse_annotation,
                                                   '_datetime': _parse_datetime,
                                                   '_list': _parse_list,
                                                   '_tile_list': _parse_tile_list})
    
    if config_file:
        if not os.path.isfile(config_file):
            raise FileNotFoundError(f"Config file {config_file} does not exist.")
    else:
        with importlib.resources.path(package='ERS_NRB.resources',
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
    processor = import_module(f'ERS_NRB.{processor_name}')
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
    for k in ['work_dir', 'scene_dir']:
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
        'annotation': 'dm,ei,id,lc,li,np,ratio',
        'logfile': 'None'
    }
    
    processing_options = {
        'acq_mode': ['IMM', 'IMP', 'APP', 'IMS', 'WSM'],
        'annotation': ['dm', 'ei', 'em', 'id', 'lc',
                       'ld', 'li', 'np', 'ratio', 'wm'],
        'dem_type': ['Copernicus 10m EEA DEM',
                     'Copernicus 30m Global DEM',
                     'Copernicus 30m Global DEM II',
                     'GETASSE30'],
        'mode': ['sar', 'nrb']}
    
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
        
        if k == 'mindate':
            v = proc_sec.get_datetime(k)
        if k == 'maxdate':
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
        if k in ['date_strict']:
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


def _parse_annotation(s):
    """Custom converter for configparser:
    https://docs.python.org/3/library/configparser.html#customizing-parser-behaviour"""
    annotation_list = _parse_list(s)
    if annotation_list is not None:
        allowed = ['dm', 'ei', 'em', 'id', 'lc', 'ld', 'li', 'np', 'ratio', 'wm']
        for layer in annotation_list:
            if layer not in allowed:
                msg = "Parameter 'annotation': Error while parsing to list; " \
                      "layer '{}' is not supported. Allowed keys:\n{}"
                raise ValueError(msg.format(layer, allowed))
    return annotation_list


def _parse_datetime(s):
    """Custom converter for configparser:
    https://docs.python.org/3/library/configparser.html#customizing-parser-behaviour"""
    if 'T' in s:
        try:
            return datetime.strptime(s, '%Y-%m-%dT%H:%M:%S')
        except ValueError as e:
            raise Exception("Parameters 'mindate/maxdate': Could not parse '{}' with datetime format "
                            "'%Y-%m-%dT%H:%M:%S'".format(s)) from e
    else:
        try:
            return datetime.strptime(s, '%Y-%m-%d')
        except ValueError as e:
            raise Exception("Parameters 'mindate/maxdate': Could not parse '{}' with datetime format "
                            "'%Y-%m-%d'".format(s)) from e


def _parse_list(s):
    """Custom converter for configparser:
    https://docs.python.org/3/library/configparser.html#customizing-parser-behaviour"""
    if s in ['', 'None']:
        return None
    else:
        return [x.strip() for x in s.split(',')]


def _parse_tile_list(s):
    """Custom converter for configparser:
    https://docs.python.org/3/library/configparser.html#customizing-parser-behaviour"""
    tile_list = _parse_list(s)
    if tile_list is not None:
        for tile in tile_list:
            if len(tile) != 5:
                raise ValueError("Parameter 'aoi_tiles': Error while parsing "
                                 "MGRS tile IDs to list; tile '{}' is not 5 "
                                 "characters long.".format(tile))
            else:
                continue
    return tile_list


def _val_cleanup(val):
    """Helper function to clean up value strings while parsing a config file."""
    return val.replace('"', '').replace("'", "")


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
