import os
import re
import time
import tempfile
from datetime import datetime, timezone
from dateutil.parser import parse as dateparse
from lxml import etree
import numpy as np
from osgeo import gdal
from spatialist import Raster, Vector, vectorize, boundary, bbox, intersect, rasterize
from spatialist.ancillary import finder
from spatialist.auxil import gdalbuildvrt, gdalwarp
from pyroSAR import identify_many, Archive
from pyroSAR.snap.util import geocode
from pyroSAR.ancillary import find_datasets
from ERS_NRB.config import get_config, geocode_conf, gdal_conf

import ERS_NRB.ancillary as ancil
import ERS_NRB.tile_extraction as tile_ex
from ERS_NRB.metadata import extract

gdal.UseExceptions()

from s1ard import dem
from s1ard.ancillary import generate_unique_id, get_max_ext
from s1ard.ard import create_rgb_vrt
from s1ard.tile_extraction import aoi_from_tile
from s1ard.metadata import stac, xml


def nrb_processing(config, product_type, scenes, datadir, outdir, tile, extent, epsg, wbm=None,
                   multithread=True, overviews=None, recursive=False, update=False):
    """
    Finalizes the generation of Sentinel-1 NRB products after the main processing steps via `pyroSAR.snap.util.geocode`
    have been executed. This includes the following:
    - Converting all GeoTIFF files to Cloud Optimized GeoTIFF (COG) format
    - Generating annotation datasets in Virtual Raster Tile (VRT) format
    - Applying the NRB product directory structure & naming convention
    - Generating metadata in XML and STAC JSON formats
    
    Parameters
    ----------
    config: dict
        Dictionary of the parsed config parameters for the current process.
    scenes: list[str]
        List of scenes to process. Either an individual scene or multiple, matching scenes (consecutive acquisitions).
    datadir: str
        The directory containing the datasets processed from the source scenes using pyroSAR.
    outdir: str
        The directory to write the final files to.
    tile: str
        ID of an MGRS tile.
    extent: dict
        Spatial extent of the MGRS tile, derived from a `spatialist.vector.Vector` object.
    epsg: int
        The CRS used for the NRB product; provided as an EPSG code.
    wbm: str, optional
        Path to a water body mask file with the dimensions of an MGRS tile.
    multithread: bool, optional
        Should `gdalwarp` use multithreading? Default is True. The number of threads used, can be adjusted in the
        config.ini file with the parameter `gdal_threads`.
    overviews: list[int], optional
        Internal overview levels to be created for each GeoTIFF file. Defaults to [2, 4, 9, 18, 36]
    recursive: bool, optional
        Find datasets recursively in the directory specified with the parameter `datadir`? Default is False.
    
    Returns
    -------
    None
    """
    # determine processing timestamp and generate unique ID
    start_time = time.time()
    proc_time = datetime.now(timezone.utc)
    t = proc_time.isoformat().encode()
    product_id = generate_unique_id(encoded_str=t)
    
    compress = config['compression']
    if overviews is None:
        overviews = [2, 4, 9, 18, 36]
    
    ids = identify_many(scenes)
    datasets = [find_datasets(directory=datadir,
                              recursive=recursive,
                              start=id.start, stop=id.start)
                for id in ids]
    if len(datasets) == 0:
        raise RuntimeError("No pyroSAR datasets were found in the directory '{}'".format(datadir))
    pattern = '[VH]{2}_gamma0-rtc'
    i = 0
    snap_dm_tile_overlap = []
    while i < len(datasets):
        pols = [x for x in datasets[i] if re.search(pattern, os.path.basename(x))]
        snap_dm_ras = re.sub(pattern, 'datamask', pols[0])
        snap_dm_vec = snap_dm_ras.replace('.tif', '.gpkg')
        
        if not all([os.path.isfile(x) for x in [snap_dm_ras, snap_dm_vec]]):
            with Raster(pols[0]) as ras:
                arr = ras.array()
                mask = ~np.isnan(arr)
                with vectorize(target=mask, reference=ras) as vec:
                    with boundary(vec, expression="value=1") as bounds:
                        if not os.path.isfile(snap_dm_ras):
                            rasterize(vectorobject=bounds, reference=ras, outname=snap_dm_ras)
                        if not os.path.isfile(snap_dm_vec):
                            bounds.write(outfile=snap_dm_vec)
        
        with Vector(snap_dm_vec) as bounds:
            with bbox(extent, epsg) as tile_geom:
                tile_geom.write(outfile=snap_dm_vec.replace('datamask', 'kml_file').replace('.gpkg', '.geojson'))
                inter = intersect(bounds, tile_geom)
                if inter is None:
                    del ids[i]
                    del datasets[i]
                else:
                    # Add snap_dm_ras to list if it overlaps with the current tile
                    snap_dm_tile_overlap.append(snap_dm_ras)
                    i += 1
                    inter.close()
    
    if len(ids) == 0:
        raise RuntimeError('None of the scenes overlap with the current tile {tile_id}: '
                           '\n{scenes}'.format(tile_id=tile, scenes=scenes))
    
    src_scenes = [i.scene for i in ids]
    
    ard_start = ids[0].start
    ard_stop = ids[0].stop
    meta = {'mission': ids[0].sensor,
            'mode': ids[0].meta['acquisition_mode'],
            'ard_spec': product_type,
            'polarization': {"['HH']": 'HH',
                             "['VV']": 'VV',
                             "['HH', 'HV']": 'HX',
                             "['HV', 'HH']": 'HX',
                             "['VV', 'HH']": 'VC',
                             "['HH', 'VV']": 'VC',
                             "['VH', 'VV']": 'VX',
                             "['VV', 'VH']": 'VX'}[str(ids[0].polarizations)],
            'start': ard_start,
            'orbitnumber': ids[0].meta['orbitNumber_abs'],
            'datatake': hex(ids[0].meta['frameNumber']).replace('x', '').upper(),
            'tile': tile,
            'id': product_id}
    
    meta_lower = dict((k, v.lower() if not isinstance(v, int) else v) for k, v in meta.items())
    skeleton_dir = '{mission}_{mode}_{ard_spec}__1S{polarization}_{start}_{orbitnumber:06}_{datatake:0>6}_{tile}_{id}'
    skeleton_files = '{mission}-{mode}-{ard_spec}-{start}-{orbitnumber:06}-{datatake:0>6}-{tile}-{suffix}.tif'
    
    ard_base = skeleton_dir.format(**meta)
    existing = finder(outdir, [ard_base.replace(product_id, '*')], foldermode=2)
    if len(existing) > 0:
        if not update:
            return 'Already processed - Skip!'
        else:
            ard_dir = existing[0]
    else:
        ard_dir = os.path.join(outdir, ard_base)
        try:
            os.makedirs(ard_dir, exist_ok=False)
        except OSError:
            return 'Already processed - Skip!'
    subdirectories = ['measurement', 'annotation', 'source', 'support']
    for subdirectory in subdirectories:
        os.makedirs(os.path.join(ard_dir, subdirectory), exist_ok=True)
    
    # nrbdir = os.path.join(outdir, skeleton_dir.format(**meta))
    # os.makedirs(nrbdir, exist_ok=True)
    
    # 'z_error': Maximum error threshold on values for LERC* compression.
    # Will be ignored if a compression algorithm is used that isn't related to LERC.
    item_map = {
        'VV_gamma0': {'suffix': 'vv-g-lin',
                      'z_error': 1e-4},
        'VH_gamma0': {'suffix': 'vh-g-lin',
                      'z_error': 1e-4},
        'HH_gamma0': {'suffix': 'hh-g-lin',
                      'z_error': 1e-4},
        'HV_gamma0': {'suffix': 'hv-g-lin',
                      'z_error': 1e-4},
        'incidenceAngleFromEllipsoid': {'suffix': 'ei',
                                        'z_error': 1e-3},
        'layoverShadowMask': {'suffix': 'dm',
                              'z_error': 0, },
        'localIncidenceAngle': {'suffix': 'li',
                                'z_error': 1e-2},
        'scatteringArea': {'suffix': 'lc',
                           'z_error': 0.1},
        'gammaSigmaRatio': {'suffix': 'gs',
                            'z_error': 1e-4},
        'acquisitionImage': {'suffix': 'id',
                             'z_error': 0},
        'VV_NESZ': {'suffix': 'np-vv',
                    'z_error': 2e-5},
        'VH_NESZ': {'suffix': 'np-vh',
                    'z_error': 2e-5},
        'HH_NESZ': {'suffix': 'np-hh',
                    'z_error': 2e-5},
        'HV_NESZ': {'suffix': 'np-hv',
                    'z_error': 2e-5}}
    
    driver = 'COG'
    ovr_resampling = 'AVERAGE'
    write_options_base = ['BLOCKSIZE=512', 'OVERVIEW_RESAMPLING={}'.format(ovr_resampling)]
    write_options = dict()
    for key in item_map:
        write_options[key] = write_options_base.copy()
        if compress is not None:
            entry = 'COMPRESS={}'.format(compress)
            write_options[key].append(entry)
            if compress.startswith('LERC'):
                entry = 'MAX_Z_ERROR={:f}'.format(item_map[key]['z_error'])
                write_options[key].append(entry)
    write_options['layoverShadowMask'].append("BIGTIFF=YES")
    write_options['acquisitionImage'].append("BIGTIFF=YES")
    
    ####################################################################################################################
    # format existing datasets found by `pyroSAR.ancillary.find_datasets`
    if len(datasets) > 1:
        files = list(zip(*datasets))
    else:
        files = datasets[0]
    
    pattern = '|'.join(item_map.keys())
    for i, item in enumerate(files):
        if isinstance(item, str):
            match = re.search(pattern, item)
            if match is not None:
                key = match.group()
            else:
                continue
        else:
            match = [re.search(pattern, x) for x in item]
            keys = [x if x is None else x.group() for x in match]
            if len(list(set(keys))) != 1:
                raise RuntimeError('file mismatch:\n{}'.format('\n'.join(item)))
            if None in keys:
                continue
            key = keys[0]
        
        if key == 'layoverShadowMask':
            # The data mask will be created later on in the processing workflow.
            continue
        
        meta_lower['suffix'] = item_map[key]['suffix']
        outname_base = skeleton_files.format(**meta_lower)
        if re.search('_gamma0', key):
            subdir = 'measurement'
        else:
            subdir = 'annotation'
        outname = os.path.join(ard_dir, subdir, outname_base)
        
        if not os.path.isfile(outname):
            os.makedirs(os.path.dirname(outname), exist_ok=True)
            bounds = [extent['xmin'], extent['ymin'], extent['xmax'], extent['ymax']]
            
            if isinstance(item, tuple):
                with Raster(list(item), list_separate=False) as ras:
                    source = ras.filename
            elif isinstance(item, str):
                source = tempfile.NamedTemporaryFile(suffix='.vrt').name
                gdalbuildvrt(item, source)
            else:
                raise TypeError('type {} is not supported: {}'.format(type(item), item))
            
            # modify temporary VRT to make sure overview levels and resampling are properly applied
            tree = etree.parse(source)
            root = tree.getroot()
            ovr = etree.SubElement(root, 'OverviewList', attrib={'resampling': ovr_resampling.lower()})
            ov = str(overviews)
            for x in ['[', ']', ',']:
                ov = ov.replace(x, '')
            ovr.text = ov
            etree.indent(root)
            tree.write(source, pretty_print=True, xml_declaration=False, encoding='utf-8')
            
            snap_nodata = 0
            gdalwarp(src=source, dst=outname, format=driver, outputBounds=bounds,
                     srcNodata=snap_nodata, dstNodata='nan', multithread=multithread,
                     creationOptions=write_options[key])
    
    if type(files[0]) == tuple:
        files = [item for tup in files for item in tup]
    gs_path = finder(ard_dir, [r'gs\.tif$'], regex=True)[0]
    measure_paths = finder(ard_dir, ['[hv]{2}-g-lin.tif$'], regex=True)
    
    # Data mask
    print("Data mask")
    if not config['dem_type'] == 'GETASSE30' and not os.path.isfile(wbm):
        raise FileNotFoundError('External water body mask could not be found: {}'.format(wbm))
    
    dm_path = gs_path.replace('-gs.tif', '-dm.tif')
    with Raster(gs_path) as ras_gs:
        extent = ras_gs.extent
    ###################################################################################################################    
    ancil.create_data_mask(outname=dm_path, valid_mask_list=snap_dm_tile_overlap, src_files=files,
                           extent=extent, epsg=epsg, driver=driver, creation_opt=write_options['layoverShadowMask'],
                           overviews=overviews, overview_resampling=ovr_resampling, wbm=wbm)
    ###################################################################################################################
    # Acquisition ID image
    with Raster(gs_path) as ras_gs:
        extent = ras_gs.extent
    ancil.create_acq_id_image(ref_tif=gs_path, valid_mask_list=snap_dm_tile_overlap, src_scenes=src_scenes,
                              extent=extent, epsg=epsg, driver=driver, creation_opt=write_options['acquisitionImage'],
                              overviews=overviews)
    ###################################################################################################################
    # sigma & gamma nought RTC
    print("sigma nought RTC")
    for item in measure_paths:
        sigma0_rtc_lin = item.replace('g-lin.tif', 's-lin.vrt')
        sigma0_rtc_log = item.replace('g-lin.tif', 's-log.vrt')
        
        if not os.path.isfile(sigma0_rtc_lin):
            ancil.vrt_pixfun(src=[item, gs_path], dst=sigma0_rtc_lin, fun='mul',
                             options={'VRTNodata': 'NaN'}, overviews=overviews, overview_resampling=ovr_resampling)
            ancil.vrt_relpath(sigma0_rtc_lin)
        
        if not os.path.isfile(sigma0_rtc_log):
            ancil.vrt_pixfun(src=sigma0_rtc_lin, dst=sigma0_rtc_log, fun='log10', scale=10,
                             options={'VRTNodata': 'NaN'}, overviews=overviews, overview_resampling=ovr_resampling)
        
        gamma0_rtc_log = item.replace('lin.tif', 'log.vrt')
        gamma0_rtc_lin = item.replace('lin.tif', 'lin.vrt')
        
        if not os.path.isfile(gamma0_rtc_lin):
            ancil.vrt_pixfun(src=[item, gs_path], dst=gamma0_rtc_lin, fun='mul',
                             options={'VRTNodata': 'NaN'}, overviews=overviews, overview_resampling=ovr_resampling)
            ancil.vrt_relpath(gamma0_rtc_lin)
        
        if not os.path.isfile(gamma0_rtc_log):
            ancil.vrt_pixfun(src=gamma0_rtc_lin, dst=gamma0_rtc_log, fun='log10', scale=10,
                             options={'VRTNodata': 'NaN'}, overviews=overviews, overview_resampling=ovr_resampling)
    ###################################################################################################################
    # log-scaled gamma nought & color composite
    if len(measure_paths) > 1:
        print("log-scaled gamma nought & color composite")
        log_vrts = []
        for item in measure_paths:
            log = item.replace('lin.tif', 'log.vrt')
            lin = item.replace('lin.tif', 'lin.vrt')
            log_vrts.append(log)
            ## The file should be already present (previous step)
            if not os.path.isfile(lin):
                ancil.vrt_pixfun(src=[item, gs_path], dst=lin, fun='mul',
                                 options={'VRTNodata': 'NaN'}, overviews=overviews, overview_resampling=ovr_resampling)
                ancil.vrt_relpath(lin)
            ## The file should be already present (previous step)
            if not os.path.isfile(log):
                ancil.vrt_pixfun(src=lin, dst=log, fun='log10', scale=10,
                                 options={'VRTNodata': 'NaN'}, overviews=overviews, overview_resampling=ovr_resampling)
        
        cc_path = re.sub('[hv]{2}', 'cc', measure_paths[0]).replace('.tif', '.vrt')
        cc_path = re.sub('[hv]{2}', 'cc', log_vrts[0])
        create_rgb_vrt(outname=cc_path, infiles=measure_paths, overviews=overviews,
                       overview_resampling=ovr_resampling)
    ####################################################################################################################
    # metadata
    print("metadata")
    nrb_tifs = finder(ard_dir, ['-[a-z]{2,3}.tif'], regex=True, recursive=True)
    meta = extract.meta_dict(config=config, target=ard_dir, src_ids=ids,
                             proc_time=proc_time, start=dateparse(ard_start),
                             stop=dateparse(ard_stop), compression=compress)
    xml.parse(meta=meta, target=ard_dir, assets=nrb_tifs, exist_ok=True)
    stac.parse(meta=meta, target=ard_dir, assets=nrb_tifs, exist_ok=True)


def main(config_file, section_name):
    update = True  # update existing products? Internal development flag.
    config = get_config(config_file=config_file, section_name=section_name)
    log = ancil.set_logging(config=config)
    geocode_prms = geocode_conf(config=config)
    gdal_prms = gdal_conf(config=config)
    
    snap_flag = True
    nrb_flag = True
    if config['mode'] == 'snap':
        nrb_flag = False
    elif config['mode'] == 'nrb':
        snap_flag = False
    
    username, password = dem.authenticate(dem_type=config['dem_type'],
                                          username=None, password=None)
    
    ####################################################################################################################
    # archive / scene selection
    scenesERS = finder(config['scene_dir'], [r'^SAR.*\.E[12]'], regex=True, recursive=True)
    scenesENVI = finder(config['scene_dir'], [r'^ASA.*\.N1'], regex=True, recursive=True)
    
    scenes = scenesERS + scenesENVI
    if not os.path.isfile(config['db_file']):
        config['db_file'] = os.path.join(config['work_dir'], config['db_file'])
    
    with Archive(dbfile=config['db_file']) as archive:
        archive.insert(scenes)
        product = 'PRI'
        if config['acq_mode'] in ['IMS']:
            product = 'SLC'
        selection = archive.select(product=product,
                                   acquisition_mode=config['acq_mode'],
                                   mindate=config['mindate'], maxdate=config['maxdate'])
    
    if len(selection) == 0:
        raise RuntimeError("No scenes could be found for acquisition mode '{acq_mode}', mindate '{mindate}' "
                           "and maxdate '{maxdate}' in directory '{scene_dir}'.".format(acq_mode=config['acq_mode'],
                                                                                        mindate=config['mindate'],
                                                                                        maxdate=config['maxdate'],
                                                                                        scene_dir=config['scene_dir']))
    log.info(f'found {len(selection)} scene(s)')
    ####################################################################################################################
    # geometry and DEM handling
    ids = identify_many(selection)
    
    geo_dict, align_dict = tile_ex.main(config=config, tr=geocode_prms['spacing'])
    # geo_dict, align_dict = tile_ex.no_aoi(ids=ids, spacing=geocode_prms['spacing'])
    aoi_tiles = list(geo_dict.keys())
    epsg_set = set([geo_dict[tile]['epsg'] for tile in list(geo_dict.keys())])
    if len(epsg_set) != 1:
        raise RuntimeError('The AOI covers multiple UTM zones: {}\n '
                           'This is currently not supported. Please refine your AOI.'.format(list(epsg_set)))
    epsg = epsg_set.pop()
    ####################################################################################################################
    # geocode & noise power - SNAP processing
    np_dict = {'sigma0': 'NESZ', 'beta0': 'NEBZ', 'gamma0': 'NEGZ'}
    np_refarea = 'sigma0'
    # for i, scene in enumerate(ids):
    #     print(scene.meta)
    if snap_flag:
        for i, scene in enumerate(ids):
            scene_base = os.path.splitext(os.path.basename(scene.scene))[0]
            out_dir_scene = os.path.join(config['out_dir'], scene_base)
            
            log.info(f'processing scene {i + 1}/{len(scenes)}: {scene.scene}')
            if os.path.isdir(out_dir_scene) and not update:
                log.info('Already processed - Skip!')
                continue
            else:
                os.makedirs(out_dir_scene, exist_ok=True)
            ########################################################################################################
            # Preparation of DEM for SAR processing
            fname_dem = os.path.join(config['tmp_dir'], scene.outname_base() + '_DEM.tif')
            if not os.path.isfile(fname_dem):
                with scene.bbox() as geom:
                    dem.mosaic(geometry=geom, outname=fname_dem,
                               dem_type=config['dem_type'],
                               username=username, password=password)
            
            log.info(f'scene {i + 1}/{len(ids)}: {scene.scene}')
            
            start_time = time.time()
            try:
                # geocoding_type = 'Range-Doppler'
                geocoding_type = 'SAR simulation cross correlation'
                geocode(infile=scene, outdir=out_dir_scene, t_srs=epsg, tmpdir=config['tmp_dir'],
                        standardGridOriginX=align_dict['xmax'], standardGridOriginY=align_dict['ymin'],
                        externalDEMFile=fname_dem, externalDEMNoDataValue=None,
                        geocoding_type=geocoding_type, **geocode_prms)
                
                t = round((time.time() - start_time), 2)
                log.info('[GEOCODE] -- {scene} -- {time}'.format(scene=scene.scene, time=t))
            except Exception as e:
                log.error('[GEOCODE] -- {scene} -- {error}'.format(scene=scene.scene, error=e))
                return
    
    ####################################################################################################################
    # NRB - final product generation
    if nrb_flag:
        # selection_grouped = groupbyTime(images=selection, function=seconds, time=60)
        selection_grouped = selection
        for t, tile in enumerate(aoi_tiles):
            outdir = os.path.join(config['out_dir'], tile)
            os.makedirs(outdir, exist_ok=True)
            
            vec = [x.geometry() for x in ids]
            extent = get_max_ext(geometries=vec)
            with bbox(coordinates=extent, crs=4326) as box:
                dem.prepare(vector=box, threads=gdal_prms['threads'],
                            dem_dir=None, wbm_dir=config['wbm_dir'],
                            dem_type=config['dem_type'],
                            tilenames=aoi_tiles, username=username, password=password,
                            dem_strict=True)
            # get the geometries of all tiles that overlap with the current scene group
            tile_vec = aoi_from_tile(tile=tile)
            del vec
            fname_wbm = os.path.join(config['wbm_dir'],
                                     config['dem_type'],
                                     '{}_WBM.tif'.format(tile_vec.mgrs))
            if not os.path.isfile(fname_wbm):
                fname_wbm = None
            
            for s, scenes in enumerate(selection_grouped):
                if isinstance(scenes, str):
                    scenes = [scenes]
                log.info('###### [NRB] Tile {t}/{t_total}: {tile} | '
                         'Scenes {s}/{s_total}: {scenes} '.format(tile=tile, t=t + 1, t_total=len(aoi_tiles),
                                                                  scenes=[os.path.basename(s) for s in scenes],
                                                                  s=s + 1, s_total=len(selection_grouped)))
                start_time = time.time()
                try:
                    nrb_processing(config=config, product_type='NRB', scenes=scenes,
                                   datadir=os.path.dirname(outdir), outdir=outdir,
                                   tile=tile, extent=geo_dict[tile]['ext'], epsg=epsg, wbm=fname_wbm,
                                   multithread=gdal_prms['multithread'], recursive=True, update=update)
                    log.info('[    NRB] -- {scenes} -- {time}'.format(scenes=scenes,
                                                                      time=round((time.time() - start_time), 2)))
                except Exception as e:
                    log.exception('[    NRB] -- {scenes} -- {error}'.format(scenes=scenes, error=e))
                    continue
        
        gdal.SetConfigOption('GDAL_NUM_THREADS', gdal_prms['threads_before'])


if __name__ == '__main__':
    main(config_file='config.ini', section_name='GENERAL')
