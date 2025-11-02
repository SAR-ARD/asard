import os
import re
import configparser
import numpy as np
import zipfile as zf
from osgeo import gdal, osr
from osgeo.gdalconst import GA_Update
from spatialist.ancillary import finder
from spatialist.auxil import gdalwarp
from spatialist.raster import Raster
from pyroSAR import identify
from pyroSAR.ancillary import Lock
from sarsen import EnvisatProduct, envisat_terrain_correction
from asard import osv
import logging

log = logging.getLogger('asard')


def get_config_keys() -> list[str]:
    return []


def get_config_section(
        parser: configparser.ConfigParser
) -> dict[str, str]:
    out = {}
    out['dem_prepare_mode'] = 'multi-UTM'
    return out


def find_datasets(
        scene: str,
        outdir: str,
        epsg: int
) -> dict[str, str] | None:
    """
    Find processed datasets for a scene in a certain CRS.

    Parameters
    ----------
    scene:
        the file name of the SAR scene
    outdir:
        the output directory in which to search for results
    epsg:
        the EPSG code defining the output projection of the processed scenes.

    Returns
    -------
        Either None if no datasets were found or a dictionary with the
        following keys and values pointing to the file names
        (polarization-specific keys depending on product availability):

         - hh-g-lin: gamma nought RTC backscatter HH polarization
         - hv-g-lin: gamma nought RTC backscatter HV polarization
         - vh-g-lin: gamma nought RTC backscatter VH polarization
         - vv-g-lin: gamma nought RTC backscatter VV polarization
         - dm: layover-shadow data mask
         - gs: gamma-sigma ratio
         - lc: local contributing area (aka scattering area)
         - li: local incident angle
    """
    basename = os.path.splitext(os.path.basename(scene))[0]
    subdir = os.path.join(outdir, basename, str(epsg))
    if not os.path.isdir(subdir):
        return
    lookup = {'dm': r'layover_shadow_mask.tif',
              'ei': r'ellipsoid_incidence_angle.tif',
              'gt': r'gt_[hv]{2}.tif',
              'gs': r'gamma_sigma_ratio.tif',
              'lc': r'lc.tif',
              'li': r'local_incidence_angle.tif'}
    out = {}
    for key, pattern in lookup.items():
        match = finder(target=subdir, matchlist=[pattern], regex=True)
        for i, item in enumerate(match):
            # add CRS if missing in metadata
            ds = gdal.Open(item, GA_Update)
            proj = ds.GetProjection()
            if proj == '':
                srs = osr.SpatialReference()
                srs.ImportFromEPSG(epsg)
                ds.SetProjection(srs.ExportToWkt())
            ds = None
            # swap North-South alignment
            with Raster(item) as ras:
                geo = ras.geo
            if geo['yres'] > 0:
                swp = item.replace('.tif', '_swp.tif')
                with Lock(swp):
                    if not os.path.isfile(swp):
                        gdalwarp(src=item, dst=swp, creationOptions='COMPRESS=ZSTD')
                match[i] = swp
        if len(match) > 0:
            if key == 'gt':
                for i, item in enumerate(match):
                    pol = re.search('[hv]{2}', item).group()
                    out[f'{pol}-g-lin'] = item
            else:
                out[key] = match[0]
    return out if len(out) > 0 else None


def lsm_encoding() -> dict[str, int]:
    """
    Get the encoding of the layover shadow mask.
    """
    return {
        'not layover, not shadow': 0,
        'layover': 1,
        'shadow': 2,
        'layover in shadow': 3,
        'nodata': 255  # dummy value
    }


def process(
        scene: str,
        outdir: str,
        dem: list[str]
) -> None:
    """
    Main function for SAR processing with sarsen.
    The processed data is in the same grid as the input DEM.
    
    Parameters
    ----------
    scene:
        The SAR scene file name.
    outdir:
        The output directory for storing the final results.
    dem:
        The DEM filenames. Can be created with :func:`cesard.dem.create`.

    Returns
    -------

    """
    basename = os.path.splitext(os.path.basename(scene))[0]
    outdir_scene = os.path.join(outdir, basename)
    os.makedirs(outdir_scene, exist_ok=True)
    
    id = identify(scene)
    osv_file = osv.get(scene=id, offline=True)
    if osv_file is not None:
        log.info(f'found OSV file: {osv_file}')
        if osv_file.endswith('.zip'):
            with zf.ZipFile(osv_file, 'r') as zip:
                zip.extractall(path=outdir_scene)
                osv_file = os.path.join(outdir_scene, zip.namelist()[0])
    else:
        raise RuntimeError(f"Could not find an OSV file for scene {os.path.basename(scene)}. "
                           f"Consider the asard.osv module for download.")
    
    for item in dem:
        with Raster(item) as ras:
            epsg = ras.epsg
        log.info(f'processing to EPSG:{epsg}')
        
        dir_out = os.path.join(outdir_scene, str(epsg))
        os.makedirs(dir_out, exist_ok=True)
        
        write_annotation = True
        for pol in id.polarizations:
            if len(id.polarizations) > 1:
                log.info(f'processing polarization {pol}')
            
            product = EnvisatProduct(
                path=scene,
                osv_file=osv_file,
                polarization='/'.join(list(pol)),
            )
            
            gt = os.path.join(dir_out, basename + f'_gt_{pol.lower()}.tif')
            
            if write_annotation:
                lc = os.path.join(dir_out, basename + '_lc.tif')
                an = dir_out
            else:
                lc = an = None
            
            rtc = envisat_terrain_correction(
                product=product,
                dem_urlpath=item,
                output_urlpath=gt,
                simulated_urlpath=lc,
                layers_urlpath=an,
                correct_radiometry='gamma_bilinear',
                interp_method='linear'
            )
            
            # modify nodata value
            ds = gdal.Open(gt, GA_Update)
            ar = ds.ReadAsArray()
            ar[~np.isfinite(ar)] = -9999
            ar[ar == 0] = -9999
            band = ds.GetRasterBand(1)
            band.SetNoDataValue(-9999)
            band.WriteArray(ar)
            band.FlushCache()
            band = None
            ds = None
            
            if lc is not None:
                ds = gdal.Open(lc, GA_Update)
                band = ds.GetRasterBand(1)
                band.SetNoDataValue(0)
                band = None
                ds = None
            
            write_annotation = False


def version_dict() -> dict[str, str]:
    import sarsen
    return {'sarsen': sarsen.__version__}
