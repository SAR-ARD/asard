import os
import re
import sys
import logging
from time import gmtime, strftime
from datetime import datetime
from lxml import etree
from textwrap import dedent
import numpy as np
from osgeo import gdal
import spatialist
from spatialist import gdalbuildvrt, Raster, bbox
import pyroSAR
from pyroSAR import identify, examine
import ERS_NRB

log = logging.getLogger('s1ard')


def vrt_pixfun(src, dst, fun, relpaths=None, scale=None, offset=None, args=None,
               options=None, overviews=None, overview_resampling=None):
    """
    Creates a VRT file for the specified source dataset(s) and adds a pixel function that should be applied on the fly
    when opening the VRT file.
    Parameters
    ----------
    src: str or list[str]
        The input dataset(s).
    dst: str
        The output dataset.
    fun: str
        A `PixelFunctionType` that should be applied on the fly when opening the VRT file. The function is applied to a
        band that derives its pixel information from the source bands. A list of possible options can be found here:
        https://gdal.org/drivers/raster/vrt.html#default-pixel-functions.
        Furthermore, the option 'decibel' can be specified, which will implement a custom pixel function that uses
        Python code for decibel conversion (10*log10).
    relpaths: bool, optional
        Should all `SourceFilename` XML elements with attribute `@relativeToVRT="0"` be updated to be paths relative to
        the output VRT file? Default is False.
    scale: int, optional
         The scale that should be applied when computing “real” pixel values from scaled pixel values on a raster band.
         Will be ignored if `fun='decibel'`.
    offset: float, optional
        The offset that should be applied when computing “real” pixel values from scaled pixel values on a raster band.
        Will be ignored if `fun='decibel'`.
    args: dict, optional
        arguments for `fun` passed as `PixelFunctionArguments`. Requires GDAL>=3.5 to be read.
    options: dict, optional
        Additional parameters passed to `gdal.BuildVRT`.
    overviews: list[int], optional
        Internal overview levels to be created for each raster file.
    overview_resampling: str, optional
        Resampling method for overview levels.
    Examples
    --------
    linear backscatter as input:
    >>> src = 's1a-iw-nrb-20220601t052704-043465-0530a1-32tpt-vh-g-lin.tif'
    decibel scaling I:
    use `log10` pixel function and additional `Scale` parameter.
    Known to display well in QGIS, but `Scale` is ignored when reading array in Python.
    >>> dst = src.replace('-lin.tif', '-log1.vrt')
    >>> create_vrt(src=src, dst=dst, fun='log10', scale=10)
    decibel scaling II:
    use custom Python pixel function. Requires additional environment variable GDAL_VRT_ENABLE_PYTHON set to YES.
    >>> dst = src.replace('-lin.tif', '-log2.vrt')
    >>> create_vrt(src=src, dst=dst, fun='decibel')
    decibel scaling III:
    use `dB` pixel function with additional `PixelFunctionArguments`. Works best but requires GDAL>=3.5.
    >>> dst = src.replace('-lin.tif', '-log3.vrt')
    >>> create_vrt(src=src, dst=dst, fun='dB', args={'fact': 10})
    """
    gdalbuildvrt(src=src, dst=dst, **options)
    tree = etree.parse(dst)
    root = tree.getroot()
    band = tree.find('VRTRasterBand')
    band.attrib['subClass'] = 'VRTDerivedRasterBand'
    
    if fun == 'decibel':
        pxfun_language = etree.SubElement(band, 'PixelFunctionLanguage')
        pxfun_language.text = 'Python'
        pxfun_type = etree.SubElement(band, 'PixelFunctionType')
        pxfun_type.text = fun
        pxfun_code = etree.SubElement(band, 'PixelFunctionCode')
        pxfun_code.text = etree.CDATA("""
    import numpy as np
    def decibel(in_ar, out_ar, xoff, yoff, xsize, ysize, raster_xsize, raster_ysize, buf_radius, gt, **kwargs):
        np.multiply(np.log10(in_ar[0], where=in_ar[0]>0.0, out=out_ar, dtype='float32'), 10.0, out=out_ar, dtype='float32')
        """)
    else:
        pixfun_type = etree.SubElement(band, 'PixelFunctionType')
        pixfun_type.text = fun
        if args is not None:
            arg = etree.SubElement(band, 'PixelFunctionArguments')
            for key, value in args.items():
                arg.attrib[key] = str(value)
        if scale is not None:
            sc = etree.SubElement(band, 'Scale')
            sc.text = str(scale)
        if offset is not None:
            off = etree.SubElement(band, 'Offset')
            off.text = str(offset)
    
    if any([overviews, overview_resampling]) is not None:
        ovr = tree.find('OverviewList')
        if ovr is None:
            ovr = etree.SubElement(root, 'OverviewList')
        if overview_resampling is not None:
            ovr.attrib['resampling'] = overview_resampling.lower()
        if overviews is not None:
            ov = str(overviews)
            for x in ['[', ']', ',']:
                ov = ov.replace(x, '')
            ovr.text = ov
    
    if relpaths:
        srcfiles = tree.xpath('//SourceFilename[@relativeToVRT="0"]')
        for srcfile in srcfiles:
            repl = os.path.relpath(srcfile.text, start=os.path.dirname(dst))
            srcfile.text = repl
            srcfile.attrib['relativeToVRT'] = '1'
    
    etree.indent(root)
    tree.write(dst, pretty_print=True, xml_declaration=False, encoding='utf-8')


def vrt_relpath(vrt):
    """
    Converts the paths to annotation datasets in a VRT file to relative paths.
    
    Parameters
    ----------
    vrt: str
        Path to the VRT file that should be updated.
    
    Returns
    -------
    None
    """
    tree = etree.parse(vrt)
    test = tree.xpath('//SourceFilename[@relativeToVRT="0"]')[0]
    repl = '../annotation/' + os.path.basename(test.text.replace('\\', '\\\\'))
    test.text = repl
    test.attrib['relativeToVRT'] = '1'
    etree.indent(tree.getroot())
    tree.write(vrt, pretty_print=True, xml_declaration=False, encoding='utf-8')


def create_data_mask(outname, valid_mask_list, src_files, extent, epsg, driver, creation_opt, overviews,
                     overview_resampling, out_format=None, wbm=None):
    """
    Creates the Data Mask file.
    
    Parameters
    ----------
    outname: str
        Full path to the output data mask file.
    valid_mask_list: list[str]
        A list of paths pointing to the datamask_ras files that intersect with the current MGRS tile.
    src_files: list[str]
        A list of paths pointing to the SNAP processed datasets of the product.
    extent: dict
        Spatial extent of the MGRS tile, derived from a `spatialist.vector.Vector` object.
    epsg: int
        The CRS used for the NRB product; provided as an EPSG code.
    driver: str
        GDAL driver to use for raster file creation.
    creation_opt: list[str]
        GDAL creation options to use for raster file creation. Should match specified GDAL driver.
    overviews: list[int]
        Internal overview levels to be created for each raster file.
    overview_resampling: str
        Resampling method for overview levels.
    out_format: str
        The desired output format of the data mask. The following options are possible:
        - 'single-layer': Individual masks will be merged into a single-layer raster file, where each mask is encoded
        as its initial bit-value.
        - 'multi-layer' (Default): Individual masks will be written into seperate bands encoded with 1 where the mask
        information is True and 0 where it is False, creating a multi-layer raster file.
    wbm: str, optional
        Path to a water body mask file with the dimensions of an MGRS tile.
    
    Returns
    -------
    None
    """
    out_nodata = 255
    
    if out_format is None:
        out_format = 'multi-layer'
    else:
        if not out_format == 'single-layer':
            raise RuntimeError("format can only be 'single-layer' or 'multi-layer'!")
    
    pols = [pol for pol in set([re.search('[VH]{2}', os.path.basename(x)).group() for x in src_files if
                                re.search('[VH]{2}', os.path.basename(x)) is not None])]
    pattern = pols[0] + '_gamma0-rtc'
    snap_gamma0 = [x for x in src_files if re.search(pattern, os.path.basename(x))]
    snap_ls_mask = [x for x in src_files if re.search('layoverShadowMask', os.path.basename(x))]
    
    dm_bands = {1: {'arr_val': 0,
                    'name': 'not layover, nor shadow'},
                2: {'arr_val': 1,
                    'name': 'layover'},
                3: {'arr_val': 2,
                    'name': 'shadow'},
                4: {'arr_val': 4,
                    'name': 'ocean water'}}
    
    tile_bounds = [extent['xmin'], extent['ymin'], extent['xmax'], extent['ymax']]
    vrt_snap_ls = '/vsimem/' + os.path.dirname(outname) + '/snap_ls.vrt'
    vrt_snap_valid = '/vsimem/' + os.path.dirname(outname) + '/snap_valid.vrt'
    vrt_snap_gamma0 = '/vsimem/' + os.path.dirname(outname) + '/snap_gamma0.vrt'
    gdalbuildvrt(snap_ls_mask, vrt_snap_ls, outputBounds=tile_bounds, void=False)
    gdalbuildvrt(valid_mask_list, vrt_snap_valid, outputBounds=tile_bounds, void=False)
    gdalbuildvrt(snap_gamma0, vrt_snap_gamma0, outputBounds=tile_bounds, void=False)
    
    with Raster(vrt_snap_ls) as ras_snap_ls:
        with bbox(extent, crs=epsg) as tile_vec:
            rows = ras_snap_ls.rows
            cols = ras_snap_ls.cols
            geotrans = ras_snap_ls.raster.GetGeoTransform()
            proj = ras_snap_ls.raster.GetProjection()
            arr_snap_dm = ras_snap_ls.array()
            
            # Add Water Body Mask
            if wbm is not None:
                with Raster(wbm) as ras_wbm:
                    
                    # Issue 1.1 Gdal somtimes creates an extra pixel, this makes sure that they are the same size
                    arr_wbm = ras_wbm.array()[:rows, :cols]
                    out_arr = np.where((arr_wbm == 1), 4, arr_snap_dm)
                    del arr_wbm
            else:
                out_arr = arr_snap_dm
                dm_bands.pop(4)
            del arr_snap_dm
            
            # Extend the shadow class of the data mask with nodata values from backscatter data and create final array
            with Raster(vrt_snap_valid)[tile_vec] as ras_snap_valid:
                with Raster(vrt_snap_gamma0)[tile_vec] as ras_snap_gamma0:
                    arr_snap_valid = ras_snap_valid.array()
                    arr_snap_gamma0 = ras_snap_gamma0.array()
                    
                    out_arr = np.nan_to_num(out_arr)
                    out_arr = np.where(((arr_snap_valid == 1) & (np.isnan(arr_snap_gamma0)) & (out_arr != 4)), 2,
                                       out_arr)
                    out_arr[np.isnan(arr_snap_valid)] = out_nodata
                    del arr_snap_gamma0
                    del arr_snap_valid
        
        if out_format == 'multi-layer':
            outname_tmp = '/vsimem/' + os.path.basename(outname) + '.vrt'
            gdriver = gdal.GetDriverByName('GTiff')
            ds_tmp = gdriver.Create(outname_tmp, cols, rows, len(dm_bands.keys()), gdal.GDT_Byte,
                                    options=['ALPHA=UNSPECIFIED', 'PHOTOMETRIC=MINISWHITE', 'BIGTIFF=YES'])
            gdriver = None
            ds_tmp.SetGeoTransform(geotrans)
            ds_tmp.SetProjection(proj)
            
            for k, v in dm_bands.items():
                band = ds_tmp.GetRasterBand(k)
                arr_val = v['arr_val']
                b_name = v['name']
                
                arr = np.full((rows, cols), 0)
                arr[out_arr == out_nodata] = out_nodata
                if arr_val == 0:
                    arr[out_arr == 0] = 1
                elif arr_val in [1, 2]:
                    arr[(out_arr == arr_val) | (out_arr == 3)] = 1
                elif arr_val == 4:
                    arr[out_arr == 4] = 1
                
                arr = arr.astype('uint8')
                band.WriteArray(arr)
                band.SetNoDataValue(out_nodata)
                band.SetDescription(b_name)
                band.FlushCache()
                band = None
                del arr
            
            ds_tmp.SetMetadataItem('TIFFTAG_DATETIME', strftime('%Y:%m:%d %H:%M:%S', gmtime()))
            ds_tmp.BuildOverviews(overview_resampling, overviews)
            outDataset_cog = gdal.GetDriverByName(driver).CreateCopy(outname, ds_tmp, strict=1, options=creation_opt)
            outDataset_cog = None
            ds_tmp = None
        elif out_format == 'single-layer':
            ras_snap_ls.write(outname, format=driver,
                              array=out_arr.astype('uint8'), nodata=out_nodata, overwrite=True, overviews=overviews,
                              options=creation_opt)
        tile_vec = None


def create_acq_id_image(ref_tif, valid_mask_list, src_scenes, extent, epsg, driver, creation_opt, overviews):
    """
    Creation of the acquisition ID image described in CARD4L NRB 2.8
    
    Parameters
    ----------
    ref_tif: str
        Full path to a reference GeoTIFF file of the product.
    valid_mask_list: list[str]
        A list of paths pointing to the datamask_ras files that intersect with the current MGRS tile.
    src_scenes: list[str]
        A list of paths pointing to the source scenes of the product.
    extent: dict
        Spatial extent of the MGRS tile, derived from a `spatialist.vector.Vector` object.
    epsg: int
        The CRS used for the NRB product; provided as an EPSG code.
    driver: str
        GDAL driver to use for raster file creation.
    creation_opt: list[str]
        GDAL creation options to use for raster file creation. Should match specified GDAL driver.
    overviews: list[int]
        Internal overview levels to be created for each raster file.
    
    Returns
    -------
    None
    """
    outname = ref_tif.replace('-gs.tif', '-id.tif')
    out_nodata = 255
    
    # If there are two source scenes, make sure that the order in the relevant lists is correct!
    if len(src_scenes) == 2:
        starts = [datetime.strptime(identify(f).start, '%Y%m%dT%H%M%S') for f in src_scenes]
        if starts[0] > starts[1]:
            src_scenes_new = [src_scenes[1]]
            src_scenes_new.append(src_scenes[0])
            src_scenes = src_scenes_new
            starts = [identify(f).start for f in src_scenes]
        start_valid = [datetime.strptime(re.search('[0-9]{8}T[0-9]{6}', os.path.basename(f)).group(),
                                         '%Y%m%dT%H%M%S') for f in valid_mask_list]
        if start_valid[0] != starts[0]:
            valid_mask_list_new = [valid_mask_list[1]]
            valid_mask_list_new.append(valid_mask_list[0])
            valid_mask_list = valid_mask_list_new
    
    tile_bounds = [extent['xmin'], extent['ymin'], extent['xmax'], extent['ymax']]
    
    arr_list = []
    for file in valid_mask_list:
        vrt_snap_valid = '/vsimem/' + os.path.dirname(outname) + 'mosaic.vrt'
        gdalbuildvrt(file, vrt_snap_valid, outputBounds=tile_bounds, void=False)
        
        with bbox(extent, crs=epsg) as tile_vec:
            with Raster(vrt_snap_valid)[tile_vec] as vrt_ras:
                vrt_arr = vrt_ras.array()
                arr_list.append(vrt_arr)
                del vrt_arr
            tile_vec = None
    
    src_scenes_clean = [os.path.basename(src).replace('.zip', '').replace('.SAFE', '') for src in src_scenes]
    tag = '{{"{src1}": 1}}'.format(src1=src_scenes_clean[0])
    out_arr = np.full(arr_list[0].shape, out_nodata)
    if len(arr_list) == 2:
        out_arr[arr_list[1] == 1] = 2
        tag = '{{"{src1}": 1, "{src2}": 2}}'.format(src1=src_scenes_clean[0], src2=src_scenes_clean[1])
    
    out_arr[arr_list[0] == 1] = 1
    creation_opt.append('TIFFTAG_IMAGEDESCRIPTION={}'.format(tag))
    
    with Raster(ref_tif) as ref_ras:
        ref_ras.write(outname, format=driver, array=out_arr.astype('uint8'), nodata=out_nodata, overwrite=True,
                      overviews=overviews, options=creation_opt)


def set_logging(config, debug=False):
    """
    Set logging for the current process.

    Parameters
    ----------
    config: dict
        Dictionary of the parsed config parameters for the current process.
    debug: bool
        Set logging level to DEBUG?

    Returns
    -------
    logging.Logger
        The log handler for the current process.
    """
    level = logging.DEBUG if debug else logging.INFO
    
    logger = logging.getLogger('ERS_NRB')
    logger.setLevel(level)
    
    log_format = "[%(asctime)s] [%(levelname)5s] %(message)s"
    formatter = logging.Formatter(fmt=log_format,
                                  datefmt='%Y-%m-%d %H:%M:%S')
    
    logfile = config['processing']['logfile']
    if logfile is not None:
        os.makedirs(os.path.dirname(logfile), exist_ok=True)
        handler = logging.FileHandler(filename=logfile, mode='a')
    else:
        handler = logging.StreamHandler(sys.stdout)
    logger.addHandler(handler)
    
    # Add header first with simple formatting
    formatter_simple = logging.Formatter("%(message)s")
    handler.setFormatter(formatter_simple)
    _log_process_config(logger=logger, config=config)
    
    # Use normal formatting from here on
    handler.setFormatter(formatter)
    
    # add pyroSAR logger
    log_pyro = logging.getLogger('pyroSAR')
    log_pyro.setLevel(level)
    log_pyro.addHandler(handler)
    
    # add s1ard logger
    log_s1ard = logging.getLogger('s1ard')
    log_s1ard.setLevel(level)
    log_s1ard.addHandler(handler)
    
    return logger


def _log_process_config(logger, config):
    """
    Adds a header to the logfile, which includes information about the current processing configuration.

    Parameters
    ----------
    logger: logging.Logger
        The logger to which the header is added to.
    config: dict
        Dictionary of the parsed config parameters for the current process.
    """
    try:
        snap_config = examine.ExamineSnap()
        core = snap_config.get_version('core')
        microwavetbx = snap_config.get_version('microwavetbx')
        snap_core = f"{core['version']} | {core['date']}"
        snap_microwavetbx = f"{microwavetbx['version']} | {microwavetbx['date']}"
    except RuntimeError:
        snap_core = 'unknown'
        snap_microwavetbx = 'unknown'
    
    header = f"""
    ====================================================================================================================
    PROCESSING CONFIGURATION

    mode                {config['processing']['mode']}
    aoi_tiles           {config['processing']['aoi_tiles']}
    aoi_geometry        {config['processing']['aoi_geometry']}
    mindate             {config['processing']['mindate'].isoformat()}
    maxdate             {config['processing']['maxdate'].isoformat()}
    acq_mode            {config['processing']['acq_mode']}

    work_dir            {config['processing']['work_dir']}
    sar_dir             {config['processing']['sar_dir']}
    tmp_dir             {config['processing']['tmp_dir']}
    wbm_dir             {config['processing']['wbm_dir']}
    
    scene_dir           {config['processing']['scene_dir']}
    db_file             {config['processing']['db_file']}
    dem_type            {config['processing']['dem_type']}
    gdal_threads        {config['processing']['gdal_threads']}

    ====================================================================================================================
    SOFTWARE

    ERS_NRB             {ERS_NRB.__version__}
    snap-core           {snap_core}
    snap-microwavetbx   {snap_microwavetbx}
    python              {sys.version}
    python-pyroSAR      {pyroSAR.__version__}
    python-spatialist   {spatialist.__version__}
    python-GDAL         {gdal.__version__}

    ====================================================================================================================
    """
    logger.info(dedent(header))
