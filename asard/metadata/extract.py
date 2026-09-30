import re
import numpy as np
from dateutil.parser import parse as dateparse
from datetime import timezone
from typing import Any
from spatialist import Raster
from spatialist.ancillary import finder
from spatialist.raster import rasterize
from pyroSAR.ancillary import get_geometry

from asard import snap
from asard.config import version_dict
from asard.metadata.mapping import ORB_MAP, NOISE_MAP, URL

from cesard.metadata.mapping import DEM_MAP, LERC_ERR_THRES
from cesard.metadata.model import (
    ARDMetadata, AxisAccuracy, BackscatterMetadata, Card4lMetadata,
    CommonMetadata, CompressionMetadata, DEMMetadata, GeometricAccuracyMetadata,
    GeometricCorrectionMetadata, GeometryMetadata, GridMetadata,
    GroundSamplingDistance, IncidenceAngleMetadata, NoiseRemovalMetadata,
    NOT_IMPLEMENTED_NUMBER, NOT_IMPLEMENTED_TEXT, PerformanceEstimate,
    ProcessingMetadata, ProductMetadata, RadiometricAccuracyMetadata,
    SourceMetadata, SourceOrbitMetadata, SourcePerformanceMetadata,
    SwathAxisMetadata,
)
from cesard.metadata.extract import (calc_enl, geometry_from_vec,
                                     calc_performance_estimates)


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
    coord_list = [
        coordinate
        for sid in src_ids
        for coordinate in sid.meta["coordinates"]
    ]
    
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
        with get_geometry(coordinates=coord_list, crs=4326) as srcvec:
            ras_srcvec = rasterize(vectorobject=srcvec, reference=ras, burn_values=[1])
            arr_srcvec = ras_srcvec.array()
            out['nodata_borderpx'] = int(np.count_nonzero(np.isnan(arr_srcvec)))
    
    proc_meta = snap.get_metadata(scene=src_ids[0], outdir=sar_dir)
    out['ML_nRgLooks'] = proc_meta['rlks'] * src_ids[0].meta['looks'][0]
    out['ML_nAzLooks'] = proc_meta['azlks'] * src_ids[0].meta['looks'][1]
    
    return out


def meta_dict(
        config: dict[str, Any],
        prod_meta: dict[str, Any],
        src_ids: list[Any],
        compression: str
) -> ARDMetadata:
    """
    Create an :class:`~cesard.metadata.model.ARDMetadata` object containing
    metadata for an ARD product and its source products.

    Parameters
    ----------
    config
        Dictionary of the parsed config parameters for the current process.
    prod_meta
        General product metadata as extracted by :func:`asard.ard.product_info`.
    src_ids
        List of :class:`~pyroSAR.drivers.ID` objects of all source products
        that overlap with the current MGRS tile.
    compression
        The compression type applied to raster files of the product.

    Returns
    -------
        The product and source metadata as an
        :class:`~cesard.metadata.model.ARDMetadata` object.
    """
    ref_tif = finder(prod_meta['dir_ard'], ['[hv]{2}-g-lin.tif$'], regex=True)[0]
    np_tifs = finder(prod_meta['dir_ard'], ['-np-[hv]{2}.tif$'], regex=True)
    prod_meta.update(get_prod_meta(
        tif=ref_tif,
        src_ids=src_ids,
        sar_dir=config['processing']['sar_dir'],
    ))
    op_mode = prod_meta['mode']
    
    src_sid = {sid.outname_base(): sid for sid in src_ids}
    sid0 = next(iter(src_sid.values()))
    
    dem_name = config['processing']['dem_type']
    dem = DEM_MAP[dem_name]
    dem_gsd = dem['gsd']
    if isinstance(dem_gsd, str):
        value, unit = re.search(
            '(^[0-9]+) (m|arcsec)$',
            dem_gsd
        ).groups()
        dem_gsd = GroundSamplingDistance(
            value=float(value),
            unit=unit,
        )
    else:
        raise ValueError('DEM GSD must be a string')
    
    # COMMON metadata
    if sid0.sensor == 'ASAR':
        constellation = 'envisat'
    elif sid0.sensor.startswith('ERS'):
        constellation = 'ers'
    else:
        raise ValueError('Unknown sensor type')
    
    pid_lookup = {'ERS1': '1', 'ERS2': '2', 'ASAR': '1'}
    psn_lookup = {'ERS1': 'ERS', 'ERS2': 'ERS', 'ASAR': 'ENVISAT'}
    platform_identifier = pid_lookup[sid0.sensor]
    platform_short_name = psn_lookup[sid0.sensor]
    platform_full_name = f'{platform_short_name.lower()}-{platform_identifier.lower()}'
    
    common = CommonMetadata(
        antenna_look_direction='RIGHT',
        constellation=constellation,
        instrument_short_name={'ERS1': 'SAR', 'ERS2': 'SAR', 'ASAR': 'ASAR'}[sid0.sensor],
        operational_mode=op_mode,
        orbit_direction={'A': 'ascending', 'D': 'descending'}[sid0.orbit],
        orbit_mean_altitude=693000.0,
        orbit_number_absolute=sid0.orbitNumber_abs,
        orbit_number_relative=sid0.orbitNumber_rel,
        platform_full_name=platform_full_name,
        platform_identifier=platform_identifier,
        platform_reference=URL['platformReference'][platform_full_name],
        platform_short_name=platform_short_name,
        polarizations=tuple(sid0.polarizations),
        polarization_mode=prod_meta['polarization'],
        processing_level='L1C',
        radar_band='C',
        radar_center_frequency=5300000000.0,
        sensor_type='RADAR',
        swath_identifier=op_mode,
        wrs_longitude_grid=sid0.meta['orbitNumber_rel'],
    )
    
    # PRODUCT metadata
    processing_center = NOT_IMPLEMENTED_TEXT
    software_versions = {
        str(key): str(value)
        for key, value in version_dict(
            processor_name=config['processing']['processor']
        ).items()
    }
    
    def _aware(value):
        if value.tzinfo is None or value.utcoffset() is None:
            return value.replace(tzinfo=timezone.utc)
        return value
    
    product = ProductMetadata(
        name='Normalised Radar Backscatter',
        product_type='NRB',
        acquisition_type='NOMINAL',
        status='PLANNED',
        access=None,
        doi=None,
        license=None,
        time_start=_aware(prod_meta['start']),
        time_stop=_aware(prod_meta['stop']),
        time_created=_aware(prod_meta['proc_time']),
        geometry=GeometryMetadata(
            bbox=tuple(float(x) for x in prod_meta['geom']['bbox']),
            bbox_native=tuple(float(x) for x in prod_meta['geom']['bbox_native']),
            geometry=dict(prod_meta['geom']['geometry']),
        ),
        grid=GridMetadata(
            epsg=prod_meta['epsg'],
            wkt=prod_meta['wkt'],
            rows=prod_meta['rows'],
            columns=prod_meta['cols'],
            pixel_spacing_row=float(prod_meta['res'][1]),
            pixel_spacing_column=float(prod_meta['res'][0]),
            transform=tuple(float(x) for x in prod_meta['transform']),
            mgrs_id=prod_meta['tile'],
            pixel_coordinate_convention='upper-left',
            number_of_border_pixels=prod_meta['nodata_borderpx'],
            definition_reference=URL['grid_definition'],
            convention_reference=URL['grid_convention'],
        ),
        processing=ProcessingMetadata(
            facility=processing_center,
            date=_aware(prod_meta['proc_time']),
            mode='PROTOTYPE',
            processor='asard',
            software=software_versions,
        ),
        card4l=Card4lMetadata(
            specification='NRB',
            version='5.5',
            document=URL['card4l_nrb'],
        ),
        compression=CompressionMetadata(
            type=compression,
            z_errors=LERC_ERR_THRES,
        ),
        dem=DEMMetadata(
            name=dem_name,
            type=dem['type'],
            reference=dem['ref'],
            access=dem['access'],
            gsd=dem_gsd,
            resampling_method='bilinear',
            egm_reference=dem['egm'],
            egm_resampling_method='bilinear' if dem['egm'] is None else None,
        ),
        backscatter=BackscatterMetadata(
            measurement='gamma0',
            convention='linear power',
            conversion_equation='10*log10(DN)',
            range_number_of_looks=float(sid0.meta['looks'][0]),
            azimuth_number_of_looks=float(sid0.meta['looks'][1]),
            equivalent_number_of_looks=calc_enl(tif=ref_tif),
        ),
        geometric_correction=GeometricCorrectionMetadata(
            algorithm=URL['geoCorrAlgorithm'],
            resampling_method='bilinear',
            accuracy=GeometricAccuracyMetadata(
                type='slant-range',
                eastern=AxisAccuracy(
                    bias=NOT_IMPLEMENTED_NUMBER,
                    standard_deviation=NOT_IMPLEMENTED_NUMBER,
                ),
                northern=AxisAccuracy(
                    bias=NOT_IMPLEMENTED_NUMBER,
                    standard_deviation=NOT_IMPLEMENTED_NUMBER,
                ),
                radial_rmse=NOT_IMPLEMENTED_NUMBER,
                reference=URL['geoCorrAccuracyReference'],
            ),
        ),
        radiometric_accuracy=RadiometricAccuracyMetadata(
            absolute=NOT_IMPLEMENTED_NUMBER,
            relative=NOT_IMPLEMENTED_NUMBER,
            reference=URL['radiometricAccuracyReference'],
        ),
        noise_removal=NoiseRemovalMetadata(
            applied=NOISE_MAP[sid0.acquisition_mode],
            algorithm=URL['noiseRemovalAlgorithm'],
        ),
        rtc_algorithm=URL['RTCAlgorithm'],
        number_of_acquisitions=len(src_sid),
        speckle_filter_applied=None,
        ellipsoidal_height=None,
        wind_normalization=None,
    )
    
    # SOURCE metadata
    if len(np_tifs) > 0:
        performance_estimates = calc_performance_estimates(files=np_tifs)
        noise_equivalent_intensity_type = 'sigma0'
    else:
        performance_estimates = {
            pol: {'minimum': None, 'mean': None, 'maximum': None}
            for pol in common.polarizations
        }
        noise_equivalent_intensity_type = None
    
    sources: dict[str, SourceMetadata] = {}
    for uid, sid in src_sid.items():
        proc_time = sid.meta['origin']['MPH']['PROC_TIME']
        proc_time = _aware(proc_time)
        try:
            proc_name, proc_version = sid.meta['origin']['MPH']['SOFTWARE_VER'].split('/')
            source_software = {proc_name: proc_version}
        except (AttributeError, ValueError):
            proc_name = None
            source_software = {}
        
        osv = sid.meta['origin']['DSD']['ORBIT STATE VECTOR 1']['FILENAME']
        sources[uid] = SourceMetadata(
            filename=sid.file,
            product_type=sid.meta['product'],
            data_geometry=sid.meta['image_geometry'].replace('_', '-').lower(),
            acquisition_type='NOMINAL',
            status='ARCHIVED',
            access=URL['source_access'],
            doi=URL['source_doi'],
            time_start=_aware(dateparse(sid.start)),
            time_stop=_aware(dateparse(sid.stop)),
            geometry=GeometryMetadata(
                bbox=tuple(float(x) for x in prod_meta['geom']['bbox']),
                geometry=dict(prod_meta['geom']['geometry']),
            ),
            processing=ProcessingMetadata(
                facility=sid.meta['origin']['MPH']['PROC_CENTER'],
                date=proc_time,
                mode='NOMINAL',
                processor=proc_name,
                software=source_software,
            ),
            orbit=SourceOrbitMetadata(
                ascending_node_date=None,
                start_time_from_ascending_node=None,
                completion_time_from_ascending_node=None,
                major_cycle_id=sid.meta['cycleNumber'],
                datatake_id=None,
                data_access=URL['orbitDataAccess'],
                data_source=ORB_MAP[sid.meta['origin']['MPH']['VECTOR_SOURCE']],
                state_vector=osv,
            ),
            swaths=(op_mode,),
            azimuth=SwathAxisMetadata(
                look_bandwidth={op_mode: NOT_IMPLEMENTED_NUMBER},
                number_of_looks={op_mode: sid.meta['looks'][1]},
                pixel_spacing={op_mode: sid.meta['spacing'][1]},
                resolution={op_mode: sid.meta['resolution'][1]},
            ),
            range=SwathAxisMetadata(
                look_bandwidth={op_mode: NOT_IMPLEMENTED_NUMBER},
                number_of_looks={op_mode: sid.meta['looks'][0]},
                pixel_spacing={op_mode: sid.meta['spacing'][0]},
                resolution={op_mode: sid.meta['resolution'][0]},
            ),
            incidence_angle=IncidenceAngleMetadata(
                minimum=float(sid.meta['incidence_nr']),
                maximum=float(sid.meta['incidence_fr']),
                mid_swath=float(sid.meta['incidence']),
            ),
            instrument_azimuth_angle=None,
            lut_applied=None,
            sensor_calibration=URL['sensorCalibration'],
            polarimetric_calibration_matrices=None,
            faraday_mean_rotation_angle=None,
            faraday_rotation_reference=URL['faradayRotationReference'],
            ionosphere_indicator=None,
            performance=SourcePerformanceMetadata(
                noise_equivalent_intensity_type=noise_equivalent_intensity_type,
                estimates={
                    pol: PerformanceEstimate(**values)
                    for pol, values in performance_estimates.items()
                },
                equivalent_number_of_looks=None,
                integrated_side_lobe_ratio=None,
                peak_side_lobe_ratio=None,
            ),
        )
    
    return ARDMetadata(
        common=common,
        product=product,
        sources=sources,
    )
