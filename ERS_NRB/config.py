import os
import configparser
from datetime import datetime
from osgeo import gdal


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
                'compression', 'db_file', 'dem_type', 'gdal_threads', 'logfile',
                'maxdate', 'mindate', 'mode', 'sar_dir', 'scene_dir', 'tmp_dir',
                'wbm_dir', 'work_dir']
    else:
        raise RuntimeError(f"unknown section: {section}. Options: 'processing'.")


def get_config(config_file):
    """Returns the content of a config file as a dictionary.
    
    Parameters
    ----------
    config_file: str
        Full path to the config file that should be parsed to a dictionary.
    
    Returns
    -------
    out_dict: dict
        Dictionary of the parsed config parameters.
    """
    
    if not os.path.isfile(config_file):
        raise FileNotFoundError("Config file {} does not exist.".format(config_file))
    
    parser = configparser.ConfigParser(allow_no_value=True,
                                       converters={'_annotation': _parse_annotation,
                                                   '_datetime': _parse_datetime,
                                                   '_list': _parse_list,
                                                   '_tile_list': _parse_tile_list})
    parser.read(config_file)
    parser_sec = parser['PROCESSING']
    
    allowed_keys = get_keys('processing')
    out_dict = {'processing': {}}
    for k, v in parser_sec.items():
        if k not in allowed_keys:
            raise ValueError("Parameter '{}' is not allowed; should be one of {}".format(k, allowed_keys))
        v = _val_cleanup(v)
        if v in ['None', 'none', '']:
            v = None
        if k == 'mode':
            allowed = ['nrb', 'sar', 'all']
            assert v in allowed, "Parameter '{}': expected to be one of {}; got '{}' instead".format(k, allowed, v)
            v = v.lower()
        if k == 'annotation':
            v = parser_sec.get_annotation(k)
        if k == 'aoi_tiles':
            if v is not None:
                v = parser_sec.get_tile_list(k)
        if k == 'aoi_geometry':
            if v is not None:
                assert os.path.isfile(v), "Parameter '{}': File {} could not be found".format(k, v)
        if k.endswith('date'):
            v = parser_sec.get_datetime(k)
        if k == 'acq_mode':
            assert v in ['IMM', 'IMP', 'APP', 'IMS', 'WSM']
        if k == 'work_dir':
            assert os.path.isdir(v), "Parameter '{}': '{}' must be an existing directory".format(k, v)
        if k.endswith('_dir') and not k == 'work_dir':
            if any(x in v for x in ['/', '\\']):
                assert os.path.isdir(v), "Parameter '{}': {} is a full path to a non-existing directory".format(k, v)
            else:
                v = os.path.join(parser_sec['work_dir'], v)
                os.makedirs(v, exist_ok=True)
        if k.endswith('_file') and not k.startswith('db'):
            if any(x in v for x in ['/', '\\']):
                assert os.path.isfile(v), "Parameter '{}': File {} could not be found".format(k, v)
            else:
                v = os.path.join(parser_sec['work_dir'], v)
                assert os.path.isfile(v), "Parameter '{}': File {} could not be found".format(k, v)
        if k == 'gdal_threads':
            v = int(v)
        if k == 'dem_type':
            allowed = ['Copernicus 10m EEA DEM', 'Copernicus 30m Global DEM II',
                       'Copernicus 30m Global DEM', 'GETASSE30', "Copernicus 90m Global DEM II"]
            assert v in allowed, "Parameter '{}': expected to be one of {}; got '{}' instead".format(k, allowed, v)
        out_dict['processing'][k] = v
    try:
        out_dict['processing']['compression']
    except:
        out_dict['processing']['compression'] = 'LZW'
    
    assert any([out_dict['processing'][k] is not None for k in ['aoi_tiles', 'aoi_geometry']])
    
    return out_dict


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


def geocode_conf(config):
    """
    Returns a dictionary of additional parameters for
    `pyroSAR.snap.util.geocode` based on processing configurations
    provided by the config file.
    
    Parameters
    ----------
    config: dict
        Dictionary of the parsed config parameters for the current process.
    
    Returns
    -------
    dict
        Dictionary of parameters that can be passed to `pyroSAR.snap.util.geocode`
    """
    return {'spacing': {'AP': 10,
                        'IM': 10,
                        'WS': 60}[config['processing']['acq_mode'][:2]],
            'export_extra': ['localIncidenceAngle', 'incidenceAngleFromEllipsoid',
                             'scatteringArea', 'layoverShadowMask', 'gammaSigmaRatio'],
            'dem_resampling_method': 'BILINEAR_INTERPOLATION',
            'img_resampling_method': 'BILINEAR_INTERPOLATION',
            'clean_edges': True,
            'clean_edges_pixels': 4,
            'cleanup': False
            }


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
