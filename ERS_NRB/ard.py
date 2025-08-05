import os
import shutil
from copy import deepcopy
from datetime import datetime, timezone
import ERS_NRB
from ERS_NRB.metadata.extract import meta_dict
from s1ard.ard import calc_product_start_stop, generate_unique_id
from s1ard.metadata import xml, stac
import logging

log = logging.getLogger('s1ard')


def append_metadata(target, config, prod_meta, src_ids, assets, compression):
    """
    Append metadata files to an ARD product.

    Parameters
    ----------
    target: str
        the path to the ARD product
    config: dict
        the configuration dictionary
    prod_meta: dict
        the product metadata as returned by :func:`product_info`
    src_ids: List[pyroSAR.drivers.ID]
        the source product objects
    assets: List[str]
        a list of assets in the ARD product
    compression: str
        the used compression algorithm

    Returns
    -------

    """
    # extract metadata
    meta = meta_dict(config=config, prod_meta=prod_meta, target=target,
                     src_ids=src_ids, compression=compression)
    
    # copy support files
    schema_dir = os.path.join(ERS_NRB.__path__[0], 'validation', 'schemas')
    if os.path.isdir(schema_dir):
        schemas = os.listdir(schema_dir)
        for schema in schemas:
            schema_in = os.path.join(schema_dir, schema)
            schema_out = os.path.join(target, 'support', schema)
            if not os.path.isfile(schema_out):
                log.info(f'creating {schema_out}')
                shutil.copy(schema_in, schema_out)
    
    # create metadata files
    xml.parse(meta=meta, target=target, assets=assets, exist_ok=True)
    stac.parse(meta=meta, target=target, assets=assets, exist_ok=True)


def product_info(product_type, src_ids, tile_id, extent, epsg):
    """
    Create ARD product metadata.

    Parameters
    ----------
    product_type: {NRB, ORB}
        the ARD product type
    src_ids: list[pyroSAR.drivers.ID]
        the source product objects
    tile_id: str
        the MGRS tile ID
    extent: dict
        the extent of the MGRS tile
    epsg: int
        the EPSG code of the MGRS tile

    Returns
    -------
    dict
        ARD product metadata
    """
    # determine processing timestamp and generate unique ID
    proc_time = datetime.now(timezone.utc)
    t = proc_time.isoformat().encode()
    product_id = generate_unique_id(encoded_str=t)
    
    sensor = src_ids[0].sensor
    acquisition_mode = src_ids[0].acquisition_mode
    if sensor in ['ERS1', 'ERS2']:
        mode = 'IM'
    elif sensor == 'ASAR':
        if acquisition_mode in ['APP', 'APS']:
            mode = 'AP'
        elif acquisition_mode in ['IMP', 'IMP']:
            mode = 'IM'
        elif acquisition_mode in ['WSM', 'WSS']:
            mode = 'WS'
        else:
            raise ValueError(f"Unknown acquisition mode: '{acquisition_mode}'")
    else:
        raise ValueError(f"Unknown sensor: '{sensor}'")
    
    ard_start, ard_stop = calc_product_start_stop(src_ids=src_ids, extent=extent, epsg=epsg)
    pol_str = '_'.join(sorted(src_ids[0].polarizations))
    meta = {'mission': sensor,
            'mode': mode,
            'phase': src_ids[0].meta['origin']['MPH']['PHASE'],
            'cycle': src_ids[0].meta['origin']['MPH']['CYCLE'],
            'product_type': product_type,
            'polarization': {'HH': 'SH',
                             'VV': 'SV',
                             'HH_HV': 'DH',
                             'VH_VV': 'DV'}[pol_str],
            'start': ard_start,
            'stop': ard_stop,
            'duration': (ard_stop - ard_start).total_seconds(),
            'proc_time': proc_time,
            'orbitnumber': src_ids[0].meta['orbitNumber_abs'],
            'orbitnumber_rel': src_ids[0].meta['orbitNumber_rel'],
            'datatake': hex(src_ids[0].meta['frameNumber']).replace('x', '').upper(),
            'tile': tile_id,
            'id': product_id}
    meta_name = deepcopy(meta)
    meta_name['mission'] = {'ERS1': 'ER1', 'ERS2': 'ER2', 'ASAR': 'ENV'}[sensor]
    meta_name['sensor'] = {'ERS1': 'S', 'ERS2': 'S', 'ASAR': 'A'}[sensor]
    meta_name['start'] = datetime.strftime(meta_name['start'], '%Y%m%dT%H%M%S')
    meta_name['duration'] = int(round(meta_name['duration']))
    del meta_name['stop']
    meta_name_lower = dict((k, v.lower() if isinstance(v, str) else v)
                           for k, v in meta_name.items())
    skeleton_dir = ('{mission}{sensor}{mode}{product_type}_{start}_{duration}__'
                    '{orbitnumber_rel:X}_S{id}_{phase}{cycle}_{polarization}')
    skeleton_files = ('{mission}{sensor}{mode}{product_type}-{start}-{duration}--'
                      '{orbitnumber_rel:x}-S{id}-{phase}{cycle}-{polarization}')
    
    meta['product_base'] = skeleton_dir.format(**meta_name)
    meta['file_base'] = skeleton_files.format(**meta_name_lower) + '-{suffix}.tif'
    return meta
