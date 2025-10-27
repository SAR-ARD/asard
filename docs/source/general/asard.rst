ARD Production
==============

The following sections give a brief overview of the major components of creating an ERS/Envisat NRB product.
All steps are comprised in function :func:`asard.processor.main`.
`asard` makes use of the `cesard <https://github.com/SAR-ARD/cesard>`_ package.

MGRS Gridding
-------------

The basis of the processing chain builds the Sentinel-2 Military Grid Reference System (MGRS) tiling system.
Hence, a reference file is needed containing the respective tile information for processing NRB products.
A KML file is available online that will be used in the following steps:

`S2A_OPER_GIP_TILPAR_MPC__20151209T095117_V20150622T000000_21000101T000000_B00.zip <https://sentiwiki.copernicus.eu/__attachments/1692737/S2A_OPER_GIP_TILPAR_MPC__20151209T095117_V20150622T000000_21000101T000000_B00.zip>`_

This file contains all relevant information about individual tiles, in particular the EPSG code of the respective UTM zone and the geometry of the tile in UTM coordinates.
This file is automatically downloaded to `~/cesard` by the function :func:`cesard.ancillary.get_kml`.
The function :func:`cesard.tile_extraction.aoi_from_tile` can be used to extract one or multiple tiles as :class:`spatialist.vector.Vector` object.

Scene Management
----------------

The L1 SAR products are managed in a local SQLite database to select scenes for processing (see pyroSAR's section on `Database Handling`_).

After loading an MGRS tile as an :class:`spatialist.vector.Vector` object and selecting all relevant overlapping scenes
from the database, processing can commence.

DEM Handling
------------

cesard offers a convenience function :func:`cesard.dem.mosaic` for creating scene-specific DEM files from various sources.
The function is based on :func:`pyroSAR.auxdata.dem_autoload` and :func:`pyroSAR.auxdata.dem_create` and will

- download all tiles of the selected source overlapping with a defined geometry
- create a GDAL VRT virtual mosaic from the tiles including gap filling over ocean areas
- create a new GeoTIFF from the VRT including geoid-ellipsoid height conversion if necessary
  (WGS84 heights are generally required for SAR processing but provided heights might be relative to a geoid like EGM2008).

SAR Processing
--------------

SNAP
====

The central function for processing backscatter data with SNAP is :func:`asard.snap.process`. It will perform all necessary steps to
generate radiometrically terrain corrected gamma/sigma naught backscatter plus all relevant additional datasets like
local incident angle and local contribution area (see argument ``export_extra``).
In a full processor run, the following functions are called in sequence:

- :func:`asard.snap.pre`: general pre-processing including

  + Orbit state vector enhancement
  + Calibration to beta naught (for RTC)

- :func:`cesard.snap.mli`: creates multi-looked image files (MLIs) per polarization if the target pixel spacing is larger than the source pixel spacing.

- :func:`cesard.snap.rtc`: radiometric terrain flattening.
  Output is backscatter in gamma naught RTC (:math:`\gamma^0_T`) and sigma naught RTC (:math:`\sigma^0_T`) as well as the scattering area (:math:`\beta^0 / \gamma^0_T`).

- :func:`cesard.snap.gsr`: computation of the gamma-sigma ratio (:math:`\sigma^0_T / \gamma^0_T`).

- :func:`cesard.snap.geo`: geocoding. This function may be called multiple times if the scene overlaps with multiple UTM zones.

The output is a BEAM-DIMAP product which consists of a `dim` metadata file and a `data` folder containing the individual image layers in ENVI format (extension `img`).
The function :func:`cesard.snap.find_datasets` can be used to collect the individual images files for a scene.

Depending on the user configuration parameters ``measurement`` and ``annotation``, some modifications to the workflow above are possible:

- :func:`cesard.snap.gsr` may be replaced by :func:`cesard.snap.sgr` to create a sigma-gamma ratio (:math:`\gamma^0_T / \sigma^0_T`)

ARD Formatting
--------------

During SAR processing, files covering a whole scene are created. In this last step, the scene-based structure is converted to the MGRS tile structure.
If one tile overlaps with multiple scenes, these scenes are first virtually mosaiced using VRT files.
The files are then subsetted to the actual tile extent, converted to Cloud Optimized GeoTIFFs (COG), and renamed to the NRB naming scheme.
All steps are performed by :func:`asard.ard.format`.
The actual file format conversion is done with :func:`spatialist.auxil.gdalwarp`, which is a simple wrapper around the gdalwarp utility of GDAL.

After all COG files have been created, GDAL VRT files are written for log scaling and conversion to other backscatter conventions using function :func:`cesard.ard.create_vrt`.

In a last step the OGC XML and STAC JSON metadata files will be written for the NRB product.

.. _Database Handling: https://pyrosar.readthedocs.io/en/latest/general/processing.html#database-handling
