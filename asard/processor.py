import os
import time
import inspect
from importlib import import_module
from osgeo import gdal
from spatialist import bbox, intersect
from spatialist.ancillary import finder
from pyroSAR import identify, identify_many, Archive

from asard.config import get_config, gdal_conf
import asard.ancillary as ancil
from asard.ard import product_info, append_metadata, get_datasets, format

gdal.UseExceptions()

from cesard import dem
from cesard.search import scene_select
from cesard.ancillary import get_max_ext, group_by_attr
import cesard.tile_extraction as tile_ex


def main(config_file, **kwargs):
    """
    Main processing function.
    
    Parameters
    ----------
    config_file: str or None
        Path to the INI configuration file.
    kwargs:
        Additional arguments overriding configuration in `config_file`.

    Returns
    -------

    """
    update = False  # update existing products? Internal development flag.
    config = get_config(config_file=config_file, **kwargs)
    log = ancil.set_logging(config=config)
    config_proc = config['processing']
    processor_name = config_proc['processor']
    processor = import_module(f'asard.{processor_name}')
    config_sar = config[processor_name]
    gdal_prms = gdal_conf(config=config)
    
    spacings = {'AP': 10, 'IM': 10, 'WS': 60}
    config_sar['spacing'] = spacings[config_proc['acq_mode'][:2]]
    
    sar_flag = 'sar' in config_proc['mode']
    nrb_flag = 'nrb' in config_proc['mode']
    
    username, password = dem.authenticate(dem_type=config_proc['dem_type'],
                                          username=None, password=None)
    
    ####################################################################################################################
    # scene selection
    log.info('collecting scenes')
    
    db_file_set = config_proc['db_file'] is not None
    scene_dir_set = config_proc['scene_dir'] is not None
    
    if db_file_set:
        archive = Archive(dbfile=config_proc['db_file'])
        if scene_dir_set:
            scenes = finder(target=config_proc['scene_dir'],
                            matchlist=[r'^SAR.*\.E[12]', r'^ASA.*\.N1'],
                            regex=True, recursive=True, foldermode=1)
            archive.insert(scenes)
    else:
        raise RuntimeError('could not select a search option. Please check your configuration.')
    
    if config_proc['scene'] is None:
        attr_search = ['sensor', 'mindate', 'maxdate',
                       'aoi_tiles', 'aoi_geometry', 'date_strict']
        dict_search = {k: config_proc[k] for k in attr_search}
        dict_search['acquisition_mode'] = config_proc['acq_mode']
        
        product = 'PRI'
        if config_proc['acq_mode'].endswith('S'):
            product = 'SLC'
        if config_proc['acq_mode'].endswith('M'):
            product = 'MR'
        dict_search['product'] = product
        
        selection, aoi_tiles = scene_select(archive=archive, **dict_search)
        
        archive.close()
        
        if len(selection) == 0:
            log.error('could not find any scenes')
            return
        
        log.info(f'found {len(selection)} scene(s)')
        scenes = identify_many(selection, sortkey='start')
    else:
        if config_proc['mode'] != ['sar']:
            raise RuntimeError("if argument 'scene' is set, the processing mode must be 'sar'")
        scenes = [identify(config_proc['scene'])]
        config_proc['acq_mode'] = scenes[0].acquisition_mode
        config_proc['product'] = scenes[0].product
        aoi_tiles = []
    
    # group scenes by absolute orbit number
    scenes_grouped = group_by_attr(scenes, lambda x: x.meta['orbitNumber_abs'])
    ####################################################################################################################
    # annotation layer selection
    annotation = config_proc['annotation']
    measurement = 'gamma'
    export_extra = None
    lookup = {'dm': 'layoverShadowMask',
              'ei': 'incidenceAngleFromEllipsoid',
              'lc': 'scatteringArea',
              'ld': 'lookDirection',
              'li': 'localIncidenceAngle',
              'np': 'NESZ',
              'gs': 'gammaSigmaRatio',
              'sg': 'sigmaGammaRatio'}
    
    if annotation is not None:
        annotation = [
            'gs' if x == 'ratio' and measurement == 'gamma'
            else 'sg' if x == 'ratio'
            else x for x in annotation]
        export_extra = []
        for layer in annotation:
            if layer in lookup:
                export_extra.append(lookup[layer])
    ####################################################################################################################
    # main SAR processing
    if sar_flag:
        for h, scenes in enumerate(scenes_grouped):
            log.info(f'SAR processing of group {h + 1}/{len(scenes_grouped)}')
            for i, scene in enumerate(scenes):
                scene_base = os.path.splitext(os.path.basename(scene.scene))[0]
                out_dir_scene = os.path.join(config_proc['sar_dir'], scene_base)
                tmp_dir_scene = os.path.join(config_proc['tmp_dir'], scene_base)
                
                log.info(f'processing scene {i + 1}/{len(scenes)}: {scene.scene}')
                if os.path.isdir(out_dir_scene) and not update:
                    log.info('Already processed - Skip!')
                    continue
                else:
                    os.makedirs(out_dir_scene, exist_ok=True)
                    os.makedirs(tmp_dir_scene, exist_ok=True)
                ########################################################################################################
                # Preparation of DEM for SAR processing
                dem_prepare_mode = config_sar['dem_prepare_mode']
                if dem_prepare_mode is not None:
                    fname_dem = dem.prepare(scene=scene, dem_type=config_proc['dem_type'],
                                            dir_out=tmp_dir_scene, username=username,
                                            password=password, mode=dem_prepare_mode,
                                            tr=(config_sar['spacing'], config_sar['spacing']))
                else:
                    fname_dem = None
                ########################################################################################################
                # main processing routine
                start_time = time.time()
                try:
                    log.info('starting SAR processing')
                    proc_args = {'scene': scene.scene,
                                 'outdir': config_proc['sar_dir'],
                                 'measurement': measurement,
                                 'tmpdir': config_proc['tmp_dir'],
                                 'dem': fname_dem,
                                 'export_extra': export_extra}
                    proc_args.update(config_sar)
                    sig = inspect.signature(processor.process)
                    accepted_params = set(sig.parameters.keys())
                    proc_args = {k: v for k, v in proc_args.items() if k in accepted_params}
                    processor.process(**proc_args)
                    t = round((time.time() - start_time), 2)
                    log.info(f'SAR processing finished in {t} seconds')
                except Exception as e:
                    log.error(msg=e)
                    raise
    ####################################################################################################################
    # NRB - final product generation
    if nrb_flag:
        product_type = 'NRB'
        log.info(f'starting {product_type} production')
        
        for s, scenes in enumerate(scenes_grouped):
            log.info(f'ARD processing of group {s + 1}/{len(scenes_grouped)}')
            log.info('preparing WBM tiles')
            vec = [x.geometry() for x in scenes]
            extent = get_max_ext(geometries=vec)
            with bbox(coordinates=extent, crs=4326) as box:
                dem.retile(vector=box, threads=gdal_prms['threads'],
                           dem_dir=None, wbm_dir=config_proc['wbm_dir'],
                           dem_type=config_proc['dem_type'],
                           tilenames=aoi_tiles, username=username, password=password,
                           dem_strict=True)
            # get the geometries of all tiles that overlap with the current scene group
            tiles = tile_ex.tile_from_aoi(vector=vec,
                                          return_geometries=True,
                                          tilenames=aoi_tiles)
            del vec
            t_total = len(tiles)
            for t, tile in enumerate(tiles):
                # select all scenes from the group whose footprint overlaps with the current tile
                scenes_sub = [x for x in scenes if intersect(tile, x.geometry())]
                scenes_sub_fnames = [x.scene for x in scenes_sub]
                outdir = os.path.join(config_proc['ard_dir'], tile.mgrs)
                os.makedirs(outdir, exist_ok=True)
                fname_wbm = os.path.join(config_proc['wbm_dir'],
                                         config_proc['dem_type'],
                                         '{}_WBM.tif'.format(tile.mgrs))
                if not os.path.isfile(fname_wbm):
                    fname_wbm = None
                add_dem = True  # add the DEM as output layer?
                dem_type = config_proc['dem_type'] if add_dem else None
                extent = tile.extent
                epsg = tile.getProjection('epsg')
                log.info(f'creating product {t + 1}/{t_total}')
                log.info(f'selected scene(s): {scenes_sub_fnames}')
                try:
                    prod_meta = product_info(product_type=product_type, src_ids=scenes_sub,
                                             tile_id=tile.mgrs, extent=extent, epsg=epsg,
                                             dir_out=outdir, update=update)
                except RuntimeError:
                    log.info('Already processed - Skip!')
                    del tiles
                    return
                log.info(f'product name: {os.path.join(outdir, prod_meta["product_base"])}')
                try:
                    src_ids, sar_assets = get_datasets(scenes=scenes_sub,
                                                       sar_dir=config_proc['sar_dir'],
                                                       extent=extent, epsg=epsg,
                                                       processor_name=config_proc['processor'])
                    
                    ard_assets = format(config=config, prod_meta=prod_meta,
                                        src_ids=src_ids, sar_assets=sar_assets,
                                        tile=tile.mgrs, extent=extent, epsg=epsg,
                                        wbm=fname_wbm, dem_type=dem_type, compress='LERC_ZSTD',
                                        multithread=gdal_prms['multithread'], annotation=annotation)
                    
                    if ard_assets is not None:
                        append_metadata(config=config, prod_meta=prod_meta,
                                        src_ids=src_ids, assets=ard_assets, compression='LERC_ZSTD')
                except Exception as e:
                    log.error(msg=e)
                    raise
            del tiles
        gdal.SetConfigOption('GDAL_NUM_THREADS', gdal_prms['threads_before'])
