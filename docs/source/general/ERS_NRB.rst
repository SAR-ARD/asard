ERS_NRB Production
==================

The following describes the current workflow for producing the ASAR/ERS Normalised Radar Backscatter product (ERS_NRB), which is being developed in study 1 of the COPA project.
This is not part of the official pyroSAR documentation.
However, as pyroSAR is the foundation of the processor, its documentation is used to outline the processor details to conveniently link to all relevant functionality.


The ASAR/ERS images are managed in a local SQLite database to select scenes for processing (see pyroSAR's section on `Database Handling`_).

After loading an MGRS tile as an :class:`spatialist.vector.Vector` object and selecting all relevant overlapping scenes
from the database, processing can commence.

The central function for processing backscatter data with SNAP is :func:`pyroSAR.snap.util.geocode`. It will perform all necessary steps to
generate radiometrically terrain corrected gamma naught backscatter plus all relevant additional datasets like
local incident angle and local contribution area (see argument ``export_extra``).

The code below presents an incomplete call to :func:`pyroSAR.snap.util.geocode` where several variables have been set implicitly.
``infile`` can either be  a single ASAR/ERS scene or multiple, which will then be mosaiced in radar geometry prior to geocoding.
``vec`` is the :class:`spatialist.vector.Vector` object
created from the S2 KML file with corner coordinate (``xmax``, ``ymin``). The resulting image tiles are aligned to this corner coordinate.

.. code-block:: python

    from pyroSAR.snap import geocode

    epsg = 32632  # just for demonstration; will be read from AOI geojson
    spacing = 10
    dem = 'Copernicus 30m Global DEM'

    geocode(infile=infile, outdir=outdir, tmpdir=tmpdir,
            t_srs=epsg, shapefile=vec, tr=spacing,
            alignToStandardGrid=True,
            standardGridOriginX=xmax, standardGridOriginY=ymin,
            demName=dem, scaling='linear',
            export_extra=['localIncidenceAngle', 'incidenceAngleFromEllipsoid',
                          'scatteringArea', 'layoverShadowMask'])


The Copernicus GLO-30 DEM can easily be used for processing as it is available via SNAP auto-download. Furthermore,
Copernicus EEA-10 DEM usage has been implemented as part of the function :func:`pyroSAR.auxdata.dem_autoload`.

Many DEMs contain heights relative to a geoid such as EGM96. For SAR processing this information needs to be converted to WGS84 ellipsoid heights.
pyroSAR offers a function :func:`pyroSAR.auxdata.get_egm_lookup` to download a conversion file used by SNAP. However, SNAP itself will also automatically download this file if not found.

Alternative to the auto-download options, a custom DEM can be passed to :func:`pyroSAR.snap.util.geocode` via argument ``externalDEMFile``.
The function :func:`pyroSAR.auxdata.dem_create` can be used to directly convert between EGM96 and WGS84 heights using GDAL.
This way, the argument ``externalDEMApplyEGM`` of function :func:`pyroSAR.snap.util.geocode` can be set to ``False`` and no additional lookup file is needed.

ASAR/ERS orbit state vector files (OSV) for enhancing the orbit location accuracy can be downloaded automatically by SNAP.
For ERS_NRB processing at least Restituted Orbit files (RESORB) are needed while the more accurate Precise Orbit Ephemerides (POEORB) delivered two weeks after scene acquisition do not provide additional benefit.

The function :func:`pyroSAR.snap.util.geocode` will create a list of plain GeoTIFF files, which are slightly larger than the actual tile to ensure full tile coverage after geocoding.
These files are then subsetted to the actual tile extent, converted to Cloud Optimized GeoTIFFs (COG), and renamed to the ERS_NRB naming scheme.
The function :func:`spatialist.auxil.gdalwarp` is used for this task, which is a simple wrapper around the gdalwarp utility of GDAL.
The following is another incomplete code example highlighting the general procedure of converting the individual images.
The ``outfile`` name is generated from information of the source images, the MGRS tile ID and the name of the respective file as written by :func:`pyroSAR.snap.util.geocode`.

.. code-block:: python

    from spatialist import gdalwarp, Raster
    from osgeo import gdal

    write_options = ['BLOCKSIZE=512',
                     'COMPRESS=LERC_ZSTD',
                     'MAX_Z_ERROR=0.001']

    with Raster(infiles, list_separate=False) as ras:
        source = ras.filename

    gdalwarp(src=source, dst=outfile,
             options={'format': 'COG',
                      'outputBounds': [xmin, ymin, xmax, ymax],
                      'creationOptions': write_options})

After all COG files have been created, GDAL VRT files are written for log scaling and sigma naught RTC backscatter computation.
The code below demonstrates the generation of a VRT file using :func:`spatialist.auxil.gdalbuildvrt` followed by an XML
modification to insert the pixel function (a way to achieve this with GDAL's gdalbuildvrt functionality has not yet been found).

.. code-block:: python

    from lxml import etree
    from spatialist import gdalbuildvrt

    def vrt_pixfun(src, dst, fun, scale=None, offset=None, options=None):
        gdalbuildvrt(src=src, dst=dst, options=options)
        tree = etree.parse(dst)
        root = tree.getroot()
        band = tree.find('VRTRasterBand')
        band.attrib['subClass'] = 'VRTDerivedRasterBand'
        pixfun = etree.SubElement(band, 'PixelFunctionType')
        pixfun.text = fun
        if scale is not None:
            sc = etree.SubElement(band, 'Scale')
            sc.text = str(scale)
        if offset is not None:
            off = etree.SubElement(band, 'Offset')
            off.text = str(offset)
        etree.indent(root)
        tree.write(dst, pretty_print=True, xml_declaration=False, encoding='utf-8')

In a last step the OGC XML and STAC JSON metadata files will be written for the ERS_NRB product.

.. _Database Handling: https://pyrosar.readthedocs.io/en/latest/general/processing.html#database-handling
