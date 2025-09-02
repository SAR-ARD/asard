import os
import re
import shutil
from spatialist.ancillary import finder
from pyroSAR import identify
from pyroSAR.snap.auxil import gpt, parse_recipe, parse_node, \
    orb_parametrize
from pyroSAR.ancillary import Lock, LockCollection

from s1ard.tile_extraction import aoi_from_scene
from s1ard.ancillary import datamask
from s1ard.config import keyval_check
from s1ard.snap import mli, rtc, gsr, sgr, geo, postprocess

# do not remove, needed for interface
from s1ard.snap import version_dict

import logging

log = logging.getLogger('s1ard')


def config_to_string(config):
    """
    Convert the values of a configuration dictionary to strings.

    Parameters
    ----------
    config: dict
        the configuration as returned by :func:`get_config_section`

    Returns
    -------
    dict
        the dictionary with the same structure but values converted to strings.
    """
    out = {}
    allowed = get_config_keys()
    for k, v in config.items():
        if k not in allowed:
            raise ValueError(f'key {k} not in allowed keys: {allowed}')
        if v is None:
            out[k] = 'None'
        elif k in ['clean_edges', 'clean_edges_pixels', 'cleanup']:
            out[k] = str(v)
        elif k == 'gpt_args' and isinstance(v, list):
            out[k] = ' '.join(v)
        else:
            out[k] = v
    return out


def get_config_keys():
    """
    Get all allowed configuration keys.

    Returns
    -------
    List[str]
    """
    return ['clean_edges', 'clean_edges_pixels', 'cleanup',
            'dem_resampling_method', 'gpt_args', 'img_resampling_method']


def get_config_section(parser, **kwargs):
    """
    Get the`config.ini` `SNAP` section content as a dictionary.

    Parameters
    ----------
    parser: configparser.ConfigParser
    kwargs: dict[str]

    Returns
    -------
    dict
    """
    out = {}
    defaults = {
        'cleanup': 'True',
        'clean_edges': 'True',
        'clean_edges_pixels': '4',
        'dem_resampling_method': 'BILINEAR_INTERPOLATION',
        'gpt_args': 'None',
        'img_resampling_method': 'BILINEAR_INTERPOLATION',
    }
    if 'SNAP' in parser.sections():
        section = parser['SNAP']
        for k, v in defaults.items():
            if k not in section and k not in kwargs:
                section[k] = v
    else:
        parser.read_dict({'SNAP': defaults})
        section = parser['SNAP']
    
    kwargs_str = config_to_string(kwargs)
    for k, v in kwargs_str.items():
        if k in get_config_keys():
            section[k] = v
    
    for k, v in section.items():
        v = keyval_check(key=k, val=v, allowed_keys=get_config_keys())
        if k == 'gpt_args':
            if v is not None:
                v = v.split(' ')
        if k == 'clean_edges_pixels':
            v = section.getint(k)
        if k in ['allow_res_osv', 'clean_edges', 'cleanup']:
            v = section.getboolean(k)
        out[k] = v
    out['dem_prepare_mode'] = 'single-4326'
    return out


def pre(src, dst, workflow, allow_res_osv=True,
        output_beta0=True, output_sigma0=True,
        output_gamma0=False, gpt_args=None):
    """
    General SAR preprocessing. The following operators are used (optional steps in brackets):
    Apply-Orbit-File->Calibration

    Parameters
    ----------
    src: str
        the file name of the source scene
    dst: str
        the file name of the target scene. Format is BEAM-DIMAP.
    workflow: str
        the output SNAP XML workflow filename.
    allow_res_osv: bool
        Also allow the less accurate RES orbit files to be used?
    output_beta0: bool
        output beta nought backscatter needed for RTC?
    output_sigma0: bool
        output sigma nought backscatter needed for NESZ?
    output_gamma0: bool
        output gamma nought backscatter needed?
    gpt_args: list[str] or None
        a list of additional arguments to be passed to the gpt call

        - e.g. ``['-x', '-c', '2048M']`` for increased tile cache size and intermediate clearing

    Returns
    -------

    See Also
    --------
    pyroSAR.snap.auxil.orb_parametrize
    """
    scene = identify(src)
    if not os.path.isfile(workflow):
        polarizations = scene.polarizations
        wf = parse_recipe('blank')
        ############################################
        read = parse_node('Read')
        read.parameters['file'] = scene.scene
        wf.insert_node(read)
        ############################################
        orb = orb_parametrize(scene=scene, formatName='ENVISAT',
                              allow_RES_OSV=allow_res_osv,
                              continueOnFail=True)
        wf.insert_node(orb, before=read.id)
        last = orb
        ############################################
        cal = parse_node('Calibration')
        wf.insert_node(cal, before=last.id)
        if scene.sensor == 'ASAR':
            if scene.acquisition_mode in ['APP', 'APS', 'IMS']:
                source_bands = [f'Intensity_{x}' for x in polarizations]
            elif scene.acquisition_mode in ['IMP', 'WSM']:
                source_bands = 'Intensity'
            else:
                raise ValueError(f'Unsupported acquisition mode: '
                                 f'{scene.acquisition_mode}')
        else:
            raise ValueError(f'Unsupported sensor: {scene.sensor}')
        cal.parameters['sourceBands'] = source_bands
        # SNAP will fail when this is set
        # cal.parameters['selectedPolarisations'] = polarizations
        cal.parameters['createBetaBand'] = False
        cal.parameters['createGammaBand'] = False
        cal.parameters['outputBetaBand'] = output_beta0
        cal.parameters['outputSigmaBand'] = output_sigma0
        cal.parameters['outputGammaBand'] = output_gamma0
        last = cal
        ############################################
        write = parse_node('Write')
        wf.insert_node(write, before=last.id)
        write.parameters['file'] = dst
        write.parameters['formatName'] = 'BEAM-DIMAP'
        ############################################
        wf.write(workflow)
    if not os.path.isfile(dst):
        gpt(xmlfile=workflow, tmpdir=os.path.dirname(dst),
            gpt_args=gpt_args)


def process(scene, outdir, measurement, spacing, dem,
            dem_resampling_method='BILINEAR_INTERPOLATION',
            img_resampling_method='BILINEAR_INTERPOLATION',
            rlks=None, azlks=None, tmpdir=None, export_extra=None,
            allow_res_osv=True, clean_edges=True, clean_edges_pixels=4,
            gpt_args=None, cleanup=True):
    """
    Main function for SAR processing with SNAP.

    Parameters
    ----------
    scene: str
        The SAR scene file name.
    outdir: str
        The output directory for storing the final results.
    measurement: {'sigma', 'gamma'}
        the backscatter measurement convention:

        - gamma: RTC gamma nought (:math:`\gamma^0_T`)
        - sigma: RTC sigma nought (:math:`\sigma^0_T`)
    spacing: int or float
        The output pixel spacing in meters.
    dem: str
        The DEM filename. Can be created with :func:`s1ard.dem.mosaic`.
    dem_resampling_method: str
        The DEM resampling method.
    img_resampling_method: str
        The image resampling method.
    rlks: int or None
        The number of range looks.
    azlks: int or None
        The number of azimuth looks.
    tmpdir: str or None
        Path to a temporary directory for intermediate products.
    export_extra: list[str] or None
        A list of ancillary layers to create. Default None: do not create any ancillary layers.
        Options:

         - DEM
         - gammaSigmaRatio: :math:`\sigma^0_T / \gamma^0_T`
         - sigmaGammaRatio: :math:`\gamma^0_T / \sigma^0_T`
         - incidenceAngleFromEllipsoid
         - layoverShadowMask
         - localIncidenceAngle
         - NESZ: noise equivalent sigma zero
         - projectedLocalIncidenceAngle
         - scatteringArea
         - lookDirection: range look direction angle
    allow_res_osv: bool
        Also allow the less accurate RES orbit files to be used?
    clean_edges: bool
        Erode noisy image edges? See :func:`pyroSAR.snap.auxil.erode_edges`.
        Does not apply to layover-shadow mask.
    clean_edges_pixels: int
        The number of pixels to erode.
    gpt_args: list[str] or None
        a list of additional arguments to be passed to the gpt call

        - e.g. ``['-x', '-c', '2048M']`` for increased tile cache size and intermediate clearing
    cleanup: bool
        Delete intermediate files after successful process termination?

    Returns
    -------

    Examples
    --------
    >>> from ERS_NRB import snap
    >>> scene = 'ASA_IMP_1PNESA20110820_092721_000000173105_00381_49533_0000.N1'
    >>> dem = 'ASA_IMP_1PNESA20110820_092721_000000173105_00381_49533_0000_DEM_GLO30.tif'
    >>> outdir = '.'
    >>> spacing = 10
    >>> export_extra = ['localIncidenceAngle', 'incidenceAngleFromEllipsoid',
    >>>                 'scatteringArea', 'layoverShadowMask', 'gammaSigmaRatio']
    >>> snap.process(scene=scene, outdir=outdir, spacing=spacing, dem=dem,
    >>>              export_extra=export_extra)
    """
    if measurement not in ['gamma', 'sigma']:
        raise RuntimeError("'measurement' must either be 'gamma' or 'sigma'")
    if export_extra is None:
        export_extra = []
    basename = os.path.splitext(os.path.basename(scene))[0]
    outdir_scene = os.path.join(outdir, basename)
    if tmpdir is None:
        tmpdir = outdir
        tmpdir_scene = os.path.join(tmpdir, basename + '_tmp')
    else:
        tmpdir_scene = os.path.join(tmpdir, basename)
    os.makedirs(outdir_scene, exist_ok=True)
    os.makedirs(tmpdir_scene, exist_ok=True)
    
    out_base = os.path.join(outdir_scene, basename)
    tmp_base = os.path.join(tmpdir_scene, basename)
    
    id = identify(scene)
    workflows = []
    
    apply_rtc = True
    ############################################################################
    # general pre-processing
    out_pre = tmp_base + '_pre.dim'
    out_pre_wf = out_pre.replace('.dim', '.xml')
    workflows.append(out_pre_wf)
    with Lock(out_pre):
        if not os.path.isfile(out_pre):
            log.info('preprocessing main scene')
            pre(src=scene, dst=out_pre, workflow=out_pre_wf,
                allow_res_osv=allow_res_osv,
                output_beta0=apply_rtc, gpt_args=gpt_args)
        else:
            log.info('main scene has already been preprocessed')
    ############################################################################
    # multi-looking
    out_mli = tmp_base + '_mli.dim'
    out_mli_wf = out_mli.replace('.dim', '.xml')
    with Lock(out_pre, soft=True):
        with Lock(out_mli):
            if not os.path.isfile(out_mli):
                mli(src=out_pre, dst=out_mli, workflow=out_mli_wf,
                    spacing=spacing, rlks=rlks, azlks=azlks, gpt_args=gpt_args)
    if not os.path.isfile(out_mli):
        out_mli = out_pre
    else:
        workflows.append(out_mli_wf)
    ############################################################################
    # radiometric terrain flattening
    out_rtc = out_gsr = out_sgr = None
    if apply_rtc:
        out_rtc = tmp_base + '_rtc.dim'
        out_rtc_wf = out_rtc.replace('.dim', '.xml')
        workflows.append(out_rtc_wf)
        output_sigma0_rtc = measurement == 'sigma' or 'gammaSigmaRatio' in export_extra
        with LockCollection([out_mli, dem], soft=True):
            with Lock(out_rtc):
                if not os.path.isfile(out_rtc):
                    log.info('radiometric terrain correction')
                    rtc(src=out_mli, dst=out_rtc, workflow=out_rtc_wf, dem=dem,
                        dem_resampling_method=dem_resampling_method,
                        sigma0=output_sigma0_rtc,
                        scattering_area='scatteringArea' in export_extra,
                        gpt_args=gpt_args)
        ########################################################################
        # gamma-sigma ratio computation
        out_gsr = None
        if 'gammaSigmaRatio' in export_extra:
            out_gsr = tmp_base + '_gsr.dim'
            out_gsr_wf = out_gsr.replace('.dim', '.xml')
            workflows.append(out_gsr_wf)
            with Lock(out_rtc, soft=True):
                with Lock(out_gsr):
                    if not os.path.isfile(out_gsr):
                        log.info('computing gamma-sigma ratio')
                        gsr(src=out_rtc, dst=out_gsr, workflow=out_gsr_wf,
                            gpt_args=gpt_args)
        ########################################################################
        # sigma-gamma ratio computation
        out_sgr = None
        if 'sigmaGammaRatio' in export_extra:
            out_sgr = tmp_base + '_sgr.dim'
            out_sgr_wf = out_sgr.replace('.dim', '.xml')
            workflows.append(out_sgr_wf)
            with Lock(out_rtc, soft=True):
                with Lock(out_sgr):
                    if not os.path.isfile(out_sgr):
                        log.info('computing sigma-gamma ratio')
                        sgr(src=out_rtc, dst=out_sgr, workflow=out_sgr_wf,
                            gpt_args=gpt_args)
    ############################################################################
    # geocoding
    
    # Process to multiple UTM zones or just one?
    # For testing purposes only.
    utm_multi = True
    
    def run():
        out_geo = out_base + '_geo_{}.dim'.format(epsg)
        out_geo_wf = out_geo.replace('.dim', '.xml')
        sources = list(filter(None, [out_mli, out_rtc, out_gsr, out_sgr]))
        with LockCollection(sources, soft=True):
            with Lock(out_geo):
                if not os.path.isfile(out_geo):
                    log.info(f'geocoding to EPSG:{epsg}')
                    scene1 = identify(out_mli)
                    pols = scene1.polarizations
                    bands0 = ['NESZ_{}'.format(pol) for pol in pols]
                    if measurement == 'gamma':
                        bands1 = ['Gamma0_{}'.format(pol) for pol in pols]
                    else:
                        bands0.extend(['Sigma0_{}'.format(pol) for pol in pols])
                        bands1 = []
                    if 'scatteringArea' in export_extra:
                        bands1.append('simulatedImage')
                    if 'lookDirection' in export_extra:
                        bands0.append('lookDirection')
                    geo(*sources,
                        dst=out_geo, workflow=out_geo_wf,
                        spacing=spacing, crs=epsg, geometry=ext,
                        export_extra=export_extra,
                        standard_grid_origin_x=align_x,
                        standard_grid_origin_y=align_y,
                        bands0=bands0, bands1=bands1, dem=dem,
                        dem_resampling_method=dem_resampling_method,
                        img_resampling_method=img_resampling_method,
                        gpt_args=gpt_args)
                    log.info('edge cleaning')
                    postprocess(out_geo, clean_edges=clean_edges,
                                clean_edges_pixels=clean_edges_pixels)
                    log.info('creating valid data masks')
                    out_geo_data = out_geo.replace('.dim', '.data')
                    pattern = r'(?:Gamma0|Sigma0)_[VH]{2}\.img$'
                    measurements = finder(out_geo_data, [pattern],
                                          regex=True, recursive=False)
                    dm_ras = os.path.join(out_geo_data, 'datamask.tif')
                    dm_vec = dm_ras.replace('.tif', '.gpkg')
                    dm_vec = datamask(measurement=measurements[0],
                                      dm_ras=dm_ras, dm_vec=dm_vec)
                else:
                    log.info(f'geocoding to EPSG:{epsg} has already been performed')
        for wf in workflows:
            wf_dst = os.path.join(outdir_scene, os.path.basename(wf))
            if wf != wf_dst and not os.path.isfile(wf_dst):
                shutil.copyfile(src=wf, dst=wf_dst)
    
    log.info('determining UTM zone overlaps')
    aois = aoi_from_scene(scene=id, multi=utm_multi)
    for aoi in aois:
        ext = aoi['extent']
        epsg = aoi['epsg']
        align_x = aoi['extent_utm']['xmin']
        align_y = aoi['extent_utm']['ymax']
        run()
    ############################################################################
    # delete intermediate files
    if cleanup:
        log.info('cleaning up')
        if id.product == 'GRD':
            # delete everything except *_pre.* products which are reused for buffering
            # this needs to be improved so that these products are also removed if they
            # are no longer needed for any buffering.
            items = finder(target=tmpdir_scene, matchlist=['*'],
                           foldermode=1, recursive=False)
            for item in items:
                if not re.search(r'_pre\.', item):
                    if os.path.isfile(item):
                        os.remove(item)
                    else:
                        shutil.rmtree(item)
        else:
            shutil.rmtree(tmpdir_scene)


def get_metadata(scene, outdir):
    """
    Get processing metadata needed for ARD product metadata.

    Parameters
    ----------
    scene: pyroSAR.drivers.ID
        the name of the SAR scene
    outdir: str
        the directory to search for processing output

    Returns
    -------
    dict
    """
    basename = scene.outname_base
    rlks = azlks = 1
    wf_mli = finder(outdir, [f'{basename}*.xml'])
    if len(wf_mli) == 1:
        wf = parse_recipe(wf_mli[0])
        if 'Multilook' in wf.operators:
            rlks = int(wf['Multilook'].parameters['nRgLooks'])
            azlks = int(wf['Multilook'].parameters['nAzLooks'])
    elif len(wf_mli) > 1:
        msg = 'found multiple multi-looking workflows:\n{}'
        raise RuntimeError(msg.format('\n'.join(wf_mli)))
    return {'azlks': azlks,
            'rlks': rlks}
