import os
import time
from osgeo import gdal
from spatialist import bbox, intersect
from spatialist.ancillary import finder
from pyroSAR import identify_many, Archive
from ERS_NRB.config import get_config, geocode_conf, gdal_conf

from ERS_NRB import snap
import ERS_NRB.ancillary as ancil
import ERS_NRB.tile_extraction as tile_ex
from ERS_NRB.ard import product_info, append_metadata

gdal.UseExceptions()

from s1ard import dem
from s1ard.ard import check_status, format, get_datasets
from s1ard.ancillary import get_max_ext, group_by_attr


def main(config_file):
    update = False  # update existing products? Internal development flag.
    config = get_config(config_file=config_file)
    log = ancil.set_logging(config=config)
    geocode_prms = geocode_conf(config=config)
    gdal_prms = gdal_conf(config=config)
    
    sar_flag = True
    nrb_flag = True
    if config['mode'] == 'sar':
        nrb_flag = False
    elif config['mode'] == 'nrb':
        sar_flag = False
    
    username, password = dem.authenticate(dem_type=config['dem_type'],
                                          username=None, password=None)
    
    ####################################################################################################################
    # archive / scene selection
    
    db_file_set = config['db_file'] is not None
    scene_dir_set = config['scene_dir'] is not None
    
    if db_file_set:
        if not os.path.isfile(config['db_file']):
            config['db_file'] = os.path.join(config['work_dir'], config['db_file'])
        archive = Archive(dbfile=config['db_file'])
        if scene_dir_set:
            scenes = finder(target=config['scene_dir'],
                            matchlist=[r'^SAR.*\.E[12]', r'^ASA.*\.N1'],
                            regex=True, recursive=True, foldermode=1)
            archive.insert(scenes)
        
        product = 'PRI'
        if config['acq_mode'].endswith('S'):
            product = 'SLC'
        selection = archive.select(product=product,
                                   acquisition_mode=config['acq_mode'],
                                   mindate=config['mindate'], maxdate=config['maxdate'])
        archive.close()
    else:
        raise RuntimeError("parameter 'db_file' is not set")
    
    if len(selection) == 0:
        msg = ("No scenes could be found for acquisition mode '{acq_mode}', "
               "mindate '{mindate}' and maxdate '{maxdate}' in directory '{scene_dir}'.")
        raise RuntimeError(msg.format(acq_mode=config['acq_mode'],
                                      mindate=config['mindate'],
                                      maxdate=config['maxdate'],
                                      scene_dir=config['scene_dir']))
    
    # group scenes by absolute orbit number
    scenes = identify_many(selection)
    scenes_grouped = group_by_attr(scenes, lambda x: x.meta['orbitNumber_abs'])
    
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
    # annotation layer selection
    annotation = config['annotation']
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
        for i, scene in enumerate(ids):
            scene_base = os.path.splitext(os.path.basename(scene.scene))[0]
            sar_dir_scene = os.path.join(config['sar_dir'], scene_base)
            
            log.info(f'processing scene {i + 1}/{len(scenes)}: {scene.scene}')
            if os.path.isdir(sar_dir_scene) and not update:
                log.info('Already processed - Skip!')
                continue
            else:
                os.makedirs(sar_dir_scene, exist_ok=True)
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
                log.info('starting SNAP processing')
                snap.process(scene=scene.scene, outdir=config['sar_dir'],
                             measurement='gamma',
                             tmpdir=config['tmp_dir'],
                             dem=fname_dem,
                             **geocode_prms)
                t = round((time.time() - start_time), 2)
                log.info(f'SNAP processing finished in {t} seconds')
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
                dem.prepare(vector=box, threads=gdal_prms['threads'],
                            dem_dir=None, wbm_dir=config['wbm_dir'],
                            dem_type=config['dem_type'],
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
                outdir = os.path.join(config['ard_dir'], tile.mgrs)
                os.makedirs(outdir, exist_ok=True)
                fname_wbm = os.path.join(config['wbm_dir'],
                                         config['dem_type'],
                                         '{}_WBM.tif'.format(tile.mgrs))
                if not os.path.isfile(fname_wbm):
                    fname_wbm = None
                add_dem = True  # add the DEM as output layer?
                dem_type = config['dem_type'] if add_dem else None
                extent = tile.extent
                epsg = tile.getProjection('epsg')
                log.info(f'creating product {t + 1}/{t_total}')
                log.info(f'selected scene(s): {scenes_sub_fnames}')
                prod_meta = product_info(product_type=product_type, src_ids=scenes_sub,
                                         tile_id=tile.mgrs, extent=extent, epsg=epsg)
                log.info(f'product name: {os.path.join(outdir, prod_meta["product_base"])}')
                
                try:
                    dir_ard = check_status(dir_out=outdir,
                                           product_base=prod_meta["product_base"],
                                           product_id=prod_meta["id"],
                                           update=update)
                except RuntimeError:
                    log.info('Already processed - Skip!')
                    del tiles
                    return
                try:
                    src_ids, sar_assets = get_datasets(scenes=scenes_sub_fnames,
                                                       sar_dir=config['sar_dir'],
                                                       extent=extent, epsg=epsg)
                    
                    ard_assets = format(config=config, prod_meta=prod_meta,
                                        src_ids=src_ids, sar_assets=sar_assets,
                                        dir_ard=dir_ard, tile=tile.mgrs, extent=extent, epsg=epsg,
                                        wbm=fname_wbm, dem_type=dem_type, compress='LERC_ZSTD',
                                        multithread=gdal_prms['multithread'], annotation=annotation)
                    
                    append_metadata(target=dir_ard, config=config, prod_meta=prod_meta,
                                    src_ids=src_ids, assets=ard_assets, compression='LERC_ZSTD')
                except Exception as e:
                    log.error(msg=e)
                    raise
            del tiles
        gdal.SetConfigOption('GDAL_NUM_THREADS', gdal_prms['threads_before'])
