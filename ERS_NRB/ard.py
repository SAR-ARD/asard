import os
import re
import shutil
from copy import deepcopy
from importlib import import_module
from datetime import datetime, timezone
from spatialist.ancillary import finder
from spatialist.raster import Raster
from spatialist.vector import Vector, intersect, bbox
from pyroSAR.drivers import ID, identify_many
from pyroSAR.ancillary import Lock
import ERS_NRB
from ERS_NRB.metadata.extract import meta_dict
from ERS_NRB.metadata.mapping import ARD_PATTERN
from s1ard.ard import calc_product_start_stop, generate_unique_id
from s1ard.metadata import xml, stac
from s1ard.ancillary import datamask
import logging

log = logging.getLogger('ERS_NRB')


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
            schema_out = os.path.join(prod_meta['dir_ard'], 'support', schema)
            if not os.path.isfile(schema_out):
                log.info(f'creating {schema_out}')
                shutil.copy(schema_in, schema_out)
    
    # create metadata files
    xml.parse(meta=meta, target=prod_meta['dir_ard'], assets=assets, exist_ok=True)
    stac.parse(meta=meta, target=prod_meta['dir_ard'], assets=assets, exist_ok=True)


def get_datasets(
        scenes: list[str],
        sar_dir: str,
        extent: dict[str, int | float],
        epsg: int,
        processor_name: str
) -> tuple[list[ID], list[dict[str, str]]]:
    """
    Collect processing output for a list of scenes. Reads metadata from all
    source products, finds matching output files in `sar_dir` and
    filters both lists depending on the actual overlap of each product's valid
    data coverage with the current MGRS tile geometry. If no output is found
    for any scene the function will raise an error.
    To obtain the extent of valid data coverage, first a binary mask raster
    file is created with the name `datamask.tif`, which is stored in the same
    folder as the processing output as found by :func:`~ERS_NRB.snap.find_datasets`.
    Then, the boundary of this binary mask is computed and stored as `datamask.gpkg`
    (see function :func:`spatialist.vector.boundary`).
    If the provided `extent` does not overlap with this boundary, the output is
    discarded. This scenario might occur when the scene's geometry read from its
    metadata overlaps with the tile but the actual extent of data does not.

    Parameters
    ----------
    scenes:
        List of scenes to process. Either an individual scene or multiple,
        matching scenes (consecutive acquisitions).
    sar_dir:
        The directory containing the SAR datasets processed from the source
        scenes using pyroSAR. The function will raise an error if the processing
        output cannot be found for all scenes in `sar_dir`.
    extent:
        Spatial extent of the MGRS tile, derived from a
        :class:`~spatialist.vector.Vector` object.
    epsg:
        The coordinate reference system as an EPSG code.
    processor_name:
        The name of the used SAR processor. The function `find_datasets` of the
        respective processor module is used.

    Returns
    -------
        List of :class:`~pyroSAR.drivers.ID` objects of all source products
        that overlap with the current MGRS tile and a list of SAR processing output
        files that match each :class:`~pyroSAR.drivers.ID` object of `ids`.
        The format of the latter is a list of dictionaries per scene with keys as
        described by e.g. :func:`ERS_NRB.snap.find_datasets`.

    See Also
    --------
    :func:`ERS_NRB.snap.find_datasets`
    """
    processor = import_module(f'ERS_NRB.{processor_name}')
    ids = identify_many(scenes, sortkey='start')
    datasets = []
    for i, _id in enumerate(ids):
        log.debug(f'collecting processing output for scene {os.path.basename(_id.scene)}')
        files = processor.find_datasets(scene=_id.scene, outdir=sar_dir, epsg=epsg)
        if files is not None:
            base = os.path.splitext(os.path.basename(_id.scene))[0]
            ocn = re.sub('(?:SLC_|GRD[FHM])_1', 'OCN__2', base)[:-5]
            # allow 1 second tolerance
            s_start = int(ocn[31])
            s_stop = int(ocn[47])
            ocn_list = list(ocn)
            s = 1
            ocn_list[31] = f'[{s_start - s}{s_start}{s_start + s}]'
            ocn_list[47] = f'[{s_stop - s}{s_stop}{s_stop + s}]'
            ocn = ''.join(ocn_list)
            log.debug(f'searching for OCN products with pattern {ocn}')
            ocn_match = finder(target=sar_dir, matchlist=[ocn], regex=True,
                               foldermode=2, recursive=False)
            if len(ocn_match) > 0:
                for v in ['owiNrcsCmod', 'owiEcmwfWindSpeed', 'owiEcmwfWindDirection']:
                    ocn_tif = os.path.join(ocn_match[0], f'{v}.tif')
                    if os.path.isfile(ocn_tif):
                        if v.endswith('Speed'):
                            files['wm_ref_speed'] = ocn_tif
                        elif v.endswith('Direction'):
                            files['wm_ref_direction'] = ocn_tif
                        else:
                            files['wm'] = ocn_tif
            datasets.append(files)
        else:
            base = os.path.basename(_id.scene)
            msg = f'cannot find processing output for scene {base} and CRS EPSG:{epsg}'
            raise RuntimeError(msg)
    
    i = 0
    while i < len(datasets):
        log.debug(f'checking tile overlap for scene {os.path.basename(ids[i].scene)}')
        measurements = [datasets[i][x] for x in datasets[i].keys()
                        if re.search('[gs]-lin', x)]
        dm_ras = os.path.join(os.path.dirname(measurements[0]), 'datamask.tif')
        dm_vec = dm_ras.replace('.tif', '.gpkg')
        dm_vec = datamask(measurement=measurements[0], dm_ras=dm_ras, dm_vec=dm_vec)
        if dm_vec is None:
            del ids[i], datasets[i]
            continue
        with Lock(dm_vec, soft=True):
            with Vector(dm_vec) as bounds:
                with bbox(extent, epsg) as tile_geom:
                    inter = intersect(bounds, tile_geom)
                    if inter is not None:
                        with Raster(dm_ras) as ras:
                            inter_min = ras.res[0] * ras.res[1]
                        if inter.getArea() < inter_min:
                            inter.close()
                            inter = None
                    if inter is None:
                        log.debug('no overlap, removing scene')
                        del ids[i]
                        del datasets[i]
                    else:
                        log.debug('overlap detected')
                        # Add dm_ras to the datasets if it overlaps with the current tile
                        datasets[i]['datamask'] = dm_ras
                        i += 1
                        inter.close()
    return ids, datasets


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
        elif acquisition_mode in ['IMP', 'IMS']:
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
                    '{orbitnumber_rel:03X}_S{id}_{phase}{cycle:03}_{polarization}')
    skeleton_files = ('{mission}{sensor}{mode}{product_type}-{start}-{duration:04}--'
                      '{orbitnumber_rel:03x}-s{id}-{phase}{cycle:03}-{polarization}')
    
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
            if existing[0] != meta['dir_ard']:
                existing_meta = re.search(ARD_PATTERN, os.path.basename(existing[0])).groupdict()
                return product_info(product_type=product_type, src_ids=src_ids,
                                    tile_id=tile_id, extent=extent, epsg=epsg,
                                    dir_out=dir_out, update=update,
                                    product_id=existing_meta['id'])
            else:
                return meta
    else:
        try:
            os.makedirs(meta['dir_ard'], exist_ok=False)
        except OSError:
            raise RuntimeError(msg)
    return meta
