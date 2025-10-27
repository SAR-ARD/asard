import os
import numpy as np
from dateutil.parser import parse as dateparse
from spatialist import Raster
from spatialist.ancillary import finder
from spatialist.raster import rasterize

import asard
from asard import snap
from asard.metadata.mapping import ORB_MAP, NOISE_MAP, URL as URL_PACKAGE

from cesard.metadata.mapping import DEM_MAP, URL as URL_BASE, LERC_ERR_THRES
from cesard.metadata.extract import calc_enl, geometry_from_vec, calc_performance_estimates, vec_from_srccoords


def get_prod_meta(tif, src_ids, sar_dir):
    """
    Collect ARD product metadata. Items are obtained from parsing the name
    of the ARD product and from reading a measurement GeoTIFF file of this
    product.
    
    Parameters
    ----------
    tif: str
        The paths to a measurement GeoTIFF file of the ARD product.
    src_ids: list[pyroSAR.drivers.ID]
        The source product object(s).
    sar_dir: str
        A path pointing to the processed SAR datasets.
    
    Returns
    -------
    dict
        A dictionary containing metadata for the product scene.
    """
    out = dict()
    coord_list = [src.meta['coordinates'] for src in src_ids]
    
    with Raster(tif) as ras:
        vec = ras.bbox()
        srs = vec.srs
        out['wkt'] = srs.ExportToWkt()
        out['epsg'] = vec.getProjection(type='epsg')
        out['rows'] = ras.rows
        out['cols'] = ras.cols
        out['res'] = ras.res
        geo = ras.geo
        out['transform'] = [geo['xres'], geo['rotation_x'], geo['xmin'],
                            geo['rotation_y'], geo['yres'], geo['ymax']]
        out['geom'] = geometry_from_vec(vectorobject=vec)
        
        # Calculate number of nodata border pixels based on source scene(s) footprint
        with vec_from_srccoords(coord_list=coord_list, crs=4326) as srcvec:
            ras_srcvec = rasterize(vectorobject=srcvec, reference=ras, burn_values=[1])
            arr_srcvec = ras_srcvec.array()
            out['nodata_borderpx'] = np.count_nonzero(np.isnan(arr_srcvec))
    
    proc_meta = snap.get_metadata(scene=src_ids[0], outdir=sar_dir)
    out['ML_nRgLooks'] = proc_meta['rlks'] * src_ids[0].meta['looks'][0]
    out['ML_nAzLooks'] = proc_meta['azlks'] * src_ids[0].meta['looks'][1]
    
    return out


def meta_dict(config, prod_meta, src_ids, compression):
    """
    Creates a dictionary containing metadata for an ARD product as well as
    its source products. The dictionary can then be used to parse XML and
    STAC JSON metadata files, respectively.
    
    Parameters
    ----------
    config: dict
        Dictionary of the parsed config parameters for the current process.
    target: str
        A path pointing to the root directory of an ARD product.
    src_ids: list[pyroSAR.drivers.ID]
        List of :class:`~pyroSAR.drivers.ID` objects of all source products
        that overlap with the current MGRS tile.
    proc_time: datetime.datetime
        The datetime object used to generate the unique product identifier.
    start: datetime.datetime
        The product acquisition start time.
    stop: datetime.datetime
        The product acquisition stop time.
    compression: str
        The compression type applied to raster files of the product.
    
    Returns
    -------
    meta: dict
        A dictionary containing an extensive collection of metadata for the
        ARD product as well as source products.
    """
    
    URL = URL_BASE
    URL.update(URL_PACKAGE)
    
    meta = {'prod': {},
            'source': {},
            'common': {}}
    
    ref_tif = finder(prod_meta['dir_ard'], ['[hv]{2}-g-lin.tif$'], regex=True)[0]
    np_tifs = finder(prod_meta['dir_ard'], ['-np-[hv]{2}.tif$'], regex=True)
    prod_meta.update(get_prod_meta(tif=ref_tif, src_ids=src_ids, sar_dir=config['processing']['sar_dir']))
    op_mode = prod_meta['mode']
    
    src_sid = {}
    for sid in src_ids:
        uid = os.path.basename(sid.scene).split('.')[0][-4:]
        src_sid[uid] = sid
    
    src0 = list(src_sid.keys())[0]  # first key/first file
    sid0 = src_sid[src0]
    # manifest0 = src_xml[src0]['manifest']
    # nsmap0 = src_xml[src0]['manifest'].nsmap
    
    dem_name = config['processing']['dem_type']
    dem_ref = DEM_MAP[dem_name]['ref']
    dem_type = DEM_MAP[dem_name]['type']
    dem_egm = DEM_MAP[dem_name]['egm']
    
    # Common metadata (sorted alphabetically)
    meta['common']['antennaLookDirection'] = 'RIGHT'
    if sid0.sensor == 'ASAR':
        constellation = 'envisat'
    elif sid0.sensor.startswith('ERS'):
        constellation = 'ers'
    else:
        raise ValueError('Unknown sensor type')
    meta['common']['antennaLookDirection'] = 'RIGHT'
    meta['common']['constellation'] = constellation
    meta['common']['instrumentShortName'] = {'ERS1': 'SAR', 'ERS2': 'SAR', 'ASAR': 'ASAR'}[sid0.sensor]
    meta['common']['operationalMode'] = op_mode
    meta['common']['orbitDirection'] = {'A': 'ascending', 'D': 'descending'}[sid0.orbit]
    meta['common']['orbitMeanAltitude'] = '{:.2e}'.format(693000)
    meta['common']['orbitNumber_abs'] = sid0.orbitNumber_abs
    meta['common']['orbitNumber_rel'] = sid0.orbitNumber_rel
    pid_lookup = {'ERS1': '1', 'ERS2': '2', 'ASAR': '1'}
    meta['common']['platformIdentifier'] = pid_lookup[sid0.sensor]
    psn_lookup = {'ERS1': 'ERS', 'ERS2': 'ERS', 'ASAR': 'ENVISAT'}
    meta['common']['platformShortName'] = psn_lookup[sid0.sensor]
    meta['common']['platformFullname'] = '{}-{}'.format(meta['common']['platformShortName'].lower(),
                                                        meta['common']['platformIdentifier'].lower())
    meta['common']['platformReference'] = URL['platformReference'][meta['common']['platformFullname']]
    meta['common']['polarisationChannels'] = sid0.polarizations
    meta['common']['polarisationMode'] = prod_meta['polarization']
    meta['common']['processingLevel'] = 'L1C'
    meta['common']['radarBand'] = 'C'
    meta['common']['radarCenterFreq'] = 5300000000
    meta['common']['sensorType'] = 'RADAR'
    meta['common']['swathIdentifier'] = op_mode
    meta['common']['wrsLongitudeGrid'] = str(sid0.meta['orbitNumber_rel'])
    
    # Product metadata (sorted alphabetically)
    meta['prod']['access'] = None
    meta['prod']['acquisitionType'] = 'NOMINAL'
    meta['prod']['ancillaryData_KML'] = URL['ancillaryData_KML']
    meta['prod']['azimuthNumberOfLooks'] = sid0.meta['looks'][1]
    meta['prod']['backscatterConvention'] = 'linear power'
    meta['prod']['backscatterConversionEq'] = '10*log10(DN)'
    meta['prod']['backscatterMeasurement'] = 'gamma0'
    meta['prod']['card4l-link'] = URL['card4l_nrb']
    meta['prod']['card4l-name'] = 'NRB'
    meta['prod']['card4l-version'] = '5.5'
    meta['prod']['compression_type'] = compression
    meta['prod']['compression_zerrors'] = LERC_ERR_THRES
    meta['prod']['crsEPSG'] = str(prod_meta['epsg'])
    meta['prod']['crsWKT'] = prod_meta['wkt']
    meta['prod']['demAccess'] = DEM_MAP[config['processing']['dem_type']]['access']
    meta['prod']['demEGMReference'] = DEM_MAP[config['processing']['dem_type']]['egm']
    meta['prod']['demEGMResamplingMethod'] = 'bilinear'
    meta['prod']['demGSD'] = DEM_MAP[config['processing']['dem_type']]['gsd']
    meta['prod']['demName'] = dem_name
    meta['prod']['demReference'] = dem_ref
    meta['prod']['demResamplingMethod'] = 'bilinear'
    meta['prod']['demType'] = dem_type
    
    meta['prod']['doi'] = None
    meta['prod']['ellipsoidalHeight'] = None
    meta['prod']['equivalentNumberLooks'] = calc_enl(tif=ref_tif)
    # meta['prod']['fileBitsPerSample'] = '32'
    # meta['prod']['fileByteOrder'] = 'little-endian'
    # meta['prod']['fileDataType'] = 'float'
    # meta['prod']['fileFormat'] = 'COG'
    # meta['prod']['filterApplied'] = False
    # meta['prod']['filterType'] = None
    # meta['prod']['filterWindowSizeCol'] = None
    # meta['prod']['filterWindowSizeLine'] = None
    
    meta['prod']['geoCorrAccuracyEasternBias'] = None
    meta['prod']['geoCorrAccuracyEasternSTDev'] = None
    meta['prod']['geoCorrAccuracyNorthernBias'] = None
    meta['prod']['geoCorrAccuracyNorthernSTDev'] = None
    meta['prod']['geoCorrAccuracy_rRMSE'] = None
    meta['prod']['geoCorrAccuracyReference'] = URL['geoCorrAccuracyReference']
    meta['prod']['geoCorrAccuracyType'] = 'slant-range'
    
    meta['prod']['geoCorrAlgorithm'] = URL['geoCorrAlgorithm']
    meta['prod']['geoCorrResamplingMethod'] = 'bilinear'
    meta['prod']['geom_stac_bbox_native'] = prod_meta['geom']['bbox_native']
    meta['prod']['geom_stac_bbox_4326'] = prod_meta['geom']['bbox']
    meta['prod']['geom_stac_geometry_4326'] = prod_meta['geom']['geometry']
    meta['prod']['geom_xml_center'] = prod_meta['geom']['center']
    meta['prod']['geom_xml_envelope'] = prod_meta['geom']['envelope']
    meta['prod']['griddingConvention'] = 'Military Grid Reference System (MGRS)'
    meta['prod']['griddingConventionURL'] = URL['griddingConventionURL']
    meta['prod']['licence'] = None
    meta['prod']['mgrsID'] = prod_meta['tile']
    meta['prod']['noiseRemovalApplied'] = NOISE_MAP[sid0.acquisition_mode]
    meta['prod']['noiseRemovalAlgorithm'] = URL['noiseRemovalAlgorithm']
    meta['prod']['numberOfAcquisitions'] = str(len(src_ids))
    meta['prod']['numBorderPixels'] = prod_meta['nodata_borderpx']
    meta['prod']['numLines'] = str(prod_meta['rows'])
    meta['prod']['numPixelsPerLine'] = str(prod_meta['cols'])
    # meta['prod']['columnSpacing'] = prod_meta['columnSpacing']
    # meta['prod']['rowSpacing'] = prod_meta['rowSpacing']
    meta['prod']['pixelCoordinateConvention'] = 'upper-left'
    meta['prod']['processingCenter'] = 'TBD'
    meta['prod']['location'] = 'TBD'
    meta['prod']['processingLevel'] = 'Level 2'
    meta['prod']['processingMode'] = 'PROTOTYPE'
    meta['prod']['processorName'] = 'asard'
    meta['prod']['processorVersion'] = asard.__version__
    meta['prod']['productName'] = 'Normalised Radar Backscatter'
    meta['prod']['productName-short'] = 'NRB'
    meta['prod']['pxSpacingColumn'] = str(prod_meta['res'][0])
    meta['prod']['pxSpacingRow'] = str(prod_meta['res'][1])
    meta['prod']['radiometricAccuracyAbsolute'] = None
    meta['prod']['radiometricAccuracyRelative'] = None
    meta['prod']['radiometricAccuracyReference'] = URL['radiometricAccuracyReference']
    meta['prod']['rangeNumberOfLooks'] = sid0.meta['looks'][0]
    meta['prod']['RTCAlgorithm'] = URL['RTCAlgorithm']
    meta['prod']['speckleFilterApplied'] = False
    meta['prod']['status'] = 'PROTOTYPE'
    meta['prod']['timeCreated'] = prod_meta['start']
    meta['prod']['timeStart'] = prod_meta['start']
    meta['prod']['timeStop'] = prod_meta['stop']
    meta['prod']['transform'] = prod_meta['transform']
    # meta['prod']['wrsLongitudeGrid'] = str(meta['common']['orbitNumbers_rel'])
    
    # Source product metadata (sorted alphabetically)
    for uid in list(src_sid.keys()):
        meta['source'][uid] = {}
        meta['source'][uid]['access'] = URL['source_access']
        meta['source'][uid]['acquisitionType'] = 'NOMINAL'
        meta['source'][uid]['ascendingNodeDate'] = None  # TBD
        # meta['source'][uid]['azimuthLookBandwidth'] = None # TBD
        meta['source'][uid]['azimuthNumberOfLooks'] = {op_mode: src_sid[uid].meta['looks'][1]}
        meta['source'][uid]['azimuthPixelSpacing'] = {op_mode: src_sid[uid].meta['spacing'][1]}
        meta['source'][uid]['azimuthResolution'] = {op_mode: src_sid[uid].meta['resolution'][1]}
        meta['source'][uid]['dataGeometry'] = src_sid[uid].meta['image_geometry'].replace('_', '-').lower()
        # meta['source'][uid]['datatakeID'] = None # needed?
        meta['source'][uid]['doi'] = URL['source_doi']
        meta['source'][uid]['faradayMeanRotationAngle'] = None
        meta['source'][uid]['faradayRotationReference'] = URL['faradayRotationReference']
        meta['source'][uid]['filename'] = src_sid[uid].file
        meta['source'][uid]['geom_stac_bbox_4326'] = prod_meta['geom']['bbox']
        meta['source'][uid]['geom_stac_geometry_4326'] = prod_meta['geom']['geometry']
        meta['source'][uid]['geom_xml_center'] = prod_meta['geom']['center']
        meta['source'][uid]['geom_xml_envelope'] = prod_meta['geom']['envelope']
        meta['source'][uid]['incidenceAngleMax'] = src_sid[uid].meta['incidence_fr']
        meta['source'][uid]['incidenceAngleMin'] = src_sid[uid].meta['incidence_nr']
        meta['source'][uid]['incidenceAngleMidSwath'] = src_sid[uid].meta['incidence']
        meta['source'][uid]['instrumentAzimuthAngle'] = None  # TBD
        meta['source'][uid]['ionosphereIndicator'] = None
        meta['source'][uid]['lutApplied'] = None  # needed?
        meta['source'][uid]['majorCycleID'] = str(sid0.meta['cycleNumber'])
        meta['source'][uid]['orbitDataAccess'] = 'https://scihub.copernicus.eu/gnss'
        meta['source'][uid]['orbitStateVector'] = src_sid[uid].meta['origin']['DSD']['ORBIT STATE VECTOR 1'][
            'FILENAME']  # Can it be more than one vector? check
        meta['source'][uid]['orbitDataSource'] = ORB_MAP[src_sid[uid].meta['origin']['MPH']['VECTOR_SOURCE']]
        if len(np_tifs) > 0:
            meta['source'][uid]['perfEstimates'] = calc_performance_estimates(files=np_tifs)
            meta['source'][uid]['perfNoiseEquivalentIntensityType'] = 'sigma0'
        else:
            stats = {stat: None for stat in ['minimum', 'mean', 'maximum']}
            pe = {pol: stats for pol in meta['common']['polarisationChannels']}
            meta['source'][uid]['perfEstimates'] = pe
            meta['source'][uid]['perfNoiseEquivalentIntensityType'] = None
        meta['source'][uid]['perfEquivalentNumberOfLooks'] = None
        meta['source'][uid]['perfIntegratedSideLobeRatio'] = None
        meta['source'][uid]['perfPeakSideLobeRatio'] = None
        meta['source'][uid]['polCalMatrices'] = None
        meta['source'][uid]['processingCenter'] = src_sid[uid].meta['origin']['MPH']['PROC_CENTER']
        proc_time = src_sid[uid].meta['origin']['MPH']['PROC_TIME']
        meta['source'][uid]['processingDate'] = proc_time
        meta['source'][uid]['processingLevel'] = 'Level 1'
        meta['source'][uid]['processingMode'] = 'NOMINAL'
        try:
            proc_name, proc_version = src_sid[uid].meta['origin']['MPH']['SOFTWARE_VER'].split('/')
            meta['source'][uid]['processorName'] = proc_name
            meta['source'][uid]['processorVersion'] = proc_version
        except:
            meta['source'][uid]['processorName'] = None
            meta['source'][uid]['processorVersion'] = None
        meta['source'][uid]['productType'] = src_sid[uid].meta['product']
        # meta['source'][uid]['rangeLookBandwidth'] = None # TBD
        meta['source'][uid]['rangeNumberOfLooks'] = {op_mode: src_sid[uid].meta['looks'][0]}
        meta['source'][uid]['rangePixelSpacing'] = {op_mode: src_sid[uid].meta['spacing'][0]}
        meta['source'][uid]['rangeResolution'] = {op_mode: src_sid[uid].meta['resolution'][0]}
        meta['source'][uid]['sensorCalibration'] = URL['sensorCalibration']
        meta['source'][uid]['status'] = 'ARCHIVED'
        meta['source'][uid]['swaths'] = [op_mode]
        meta['source'][uid]['timeCompletionFromAscendingNode'] = None
        meta['source'][uid]['timeStartFromAscendingNode'] = None
        meta['source'][uid]['timeStart'] = dateparse(src_sid[uid].start)
        meta['source'][uid]['timeStop'] = dateparse(src_sid[uid].stop)
    return meta
