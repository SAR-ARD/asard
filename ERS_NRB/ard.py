import os
import re
import shutil
from copy import deepcopy
from datetime import datetime, timezone
from spatialist.ancillary import finder
import ERS_NRB
from ERS_NRB.metadata.extract import meta_dict
from ERS_NRB.metadata.mapping import ARD_PATTERN
from s1ard.ard import calc_product_start_stop, generate_unique_id
from s1ard.metadata import xml, stac
import logging

log = logging.getLogger('s1ard')


def append_metadata(config, prod_meta, src_ids, assets, compression):
    """
    Append metadata files to an ARD product.

    Parameters
    ----------
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
    meta = meta_dict(config=config, prod_meta=prod_meta,
                     src_ids=src_ids, compression=compression)
    
    # copy support files
    schema_dir = os.path.join(ERS_NRB.__path__[0], 'validation', 'schemas')
    if os.path.isdir(schema_dir):
        schemas = os.listdir(schema_dir)
        for schema in schemas:
            schema_in = os.path.join(schema_dir, schema)
            schema_out = os.path.join(prod_meta['ard_dir'], 'support', schema)
            if not os.path.isfile(schema_out):
                log.info(f'creating {schema_out}')
                shutil.copy(schema_in, schema_out)
    
    # create metadata files
    xml.parse(meta=meta, target=prod_meta['ard_dir'], assets=assets, exist_ok=True)
    stac.parse(meta=meta, target=prod_meta['ard_dir'], assets=assets, exist_ok=True)


def product_info(product_type, src_ids, tile_id, extent, epsg,
                 dir_out, update=False, product_id=None):
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
    if product_id is None:
        t = proc_time.isoformat().encode()
        product_id = generate_unique_id(encoded_str=t, length=3)
    
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
    skeleton_dir = ('{mission}{sensor}{mode}{product_type}_{start}_{duration:04}__'
                    '{orbitnumber_rel:X}_S{id}_{phase}{cycle}_{polarization}')
    skeleton_files = ('{mission}{sensor}{mode}{product_type}-{start}-{duration:04}--'
                      '{orbitnumber_rel:x}-S{id}-{phase}{cycle}-{polarization}')
    
    meta['product_base'] = skeleton_dir.format(**meta_name)
    meta['dir_ard'] = os.path.join(dir_out, meta['product_base'])
    meta['file_base'] = skeleton_files.format(**meta_name_lower) + '-{suffix}.tif'
    
    # check existence of products
    msg = 'Already processed - Skip!'
    pattern = meta['product_base'].replace(product_id, '*')
    existing = finder(dir_out, [pattern], foldermode=2)
    if len(existing) > 0:
        if not update:
            raise RuntimeError(msg)
        else:
            existing_meta = re.search(ARD_PATTERN, os.path.basename(existing[0])).groupdict()
            product_info(product_type=product_type, src_ids=src_ids,
                         tile_id=tile_id, extent=extent, epsg=epsg,
                         dir_out=dir_out, update=update,
                         product_id=existing_meta['id'])
    else:
        try:
            os.makedirs(meta['dir_ard'], exist_ok=False)
        except OSError:
            raise RuntimeError(msg)
    return meta
