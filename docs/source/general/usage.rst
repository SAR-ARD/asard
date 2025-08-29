Usage
=====

Configuration
-------------
Usage of the ERS_NRB package currently relies on a configuration file that needs to be set up by the user. The configuration
file follows the INI format, which uses plain text to store properties as key-value pairs. INI files can be created and
opened with any text editor. An example ``config.ini`` file for the ERS_NRB package can be found here:

https://github.com/SAR-ARD/ERS_NRB/blob/main/config.ini

The following provides an overview of the parameters the ``config.ini`` should contain and anything that should be
considered when selecting their values:

- **mode**

Options: ``all | nrb | snap``

This parameter determines if the entire processing chain should be executed or only part of it.

- **aoi_tiles** & **aoi_geometry**

The area of interest (AOI) for which ERS-NRB products should be created.

``aoi_tiles`` can be used to define the area of interest via MGRS tile IDs, which must be provided comma-separated (e.g.,
``aoi_tiles = 32TNS, 32TMT, 32TMS``). ``aoi_geometry`` defines the area of interest via a full path to a vector file
supported by :class:`spatialist.vector.Vector`. This option will automatically search for overlapping MGRS tiles and use
these for processing.
Both parameters are optional and can be set to ``None`` or left empty. ``aoi_tiles`` overrides ``aoi_geometry``.
If neither is defined, all tiles overlapping with the scene search result are processed.

- **mindate** & **maxdate**

The time period to create ERS-NRB products for. Allowed date formats are ``%Y-%m-%d`` and ``%Y-%m-%dT%H:%M:%S``.

- **acq_mode**

Options: ``IMP | IMM | IMS | APP | WSM``

The acquisition mode of the source scenes that should be processed.

- **work_dir** & **scene_dir**

Both need to be provided as full paths to existing directories. ``work_dir`` is the main directory in which any
subdirectories and files are stored that are generated during processing. ``scene_dir`` will be searched recursively for
any ASAR/ERS scenes using the regex pattern ``'^SAR.*\.E[12]$'`` and ``'^ASA.*\.N1$'``.

- **tmp_dir**, **dem_dir** & **wbm_dir**

Processing ERS-NRB products creates many intermediate files that are expected to be stored in separate subdirectories. The
default values provided in the example configuration file linked above are recommended and will automatically create
subdirectories relative to the directory specified with ``work_dir``. E.g., ``dem_dir = DEM`` will create the subdirectory
``/<work_dir>/DEM``. Optionally, full paths to existing directories can be provided for all of these parameters.

- **db_file**

Any ASAR/ERS scenes found in ``scene_dir`` will be stored in a database file created by :class:`pyrosar.drivers.Archive`.
With ``db_file`` either a full path to an existing database can be provided or it will be created in ``work_dir`` if only
a filename is provided. E.g., ``db_file = scenes.db`` will automatically create the database file ``/<work_dir>/scenes.db``.

- **kml_file**

The Sentinel-2 Military Grid Reference System (MGRS) tiling system establishes the basis of the processing chain and a
local reference file containing the respective tile information for processing S1-NRB products is needed. The official
KML file provided by ESA can be retrieved `here <https://sentinel.esa.int/documents/247904/1955685/S2A_OPER_GIP_TILPAR_MPC__20151209T095117_V20150622T000000_21000101T000000_B00.kml>`_.
With the ``kml_file`` parameter either a full path to this reference file can be provided or it is expected to be located
in the directory provided with ``work_dir`` if only a filename is provided. E.g., the processor expects to find
``/<work_dir>/s2_grid.kml`` if ``kml_file = s2_grid.kml``.

- **dem_type**

Options: ``Copernicus 10m EEA DEM | Copernicus 30m Global DEM II | Copernicus 30m Global DEM | GETASSE30``

The Digital Elevation Model (DEM) that should be used for processing.

Note that water body masks are not available for "Copernicus 30m Global DEM" and "GETASSE30", and will therefore not be
included in the product data mask. "Copernicus 10m EEA DEM" and "Copernicus 30m Global DEM II" (both include water body masks)
are retrieved from the `Copernicus Space Component Data Access system (CSCDA) <https://spacedata.copernicus.eu/web/cscda/data-access/registration>`_,
which requires registration. You will need to set up the environment variables ``FTP_USER`` and ``FTP_PASS`` to authenticate.

- **gdal_threads**

Temporarily changes GDAL_NUM_THREADS during processing. Will be reset after processing has finished.

- **compression**

Selects the compression of the final product. LZW if not specified. See https://gdal.org/drivers/raster/gtiff.html#creation-options for options.

Sections
^^^^^^^^
Configuration files in INI format can have different sections. Each section begins at a section name and ends at the next
section name. The ``config.ini`` file used with the ERS_NRB package should at least have a dedicated section for processing
related parameters. This section is by default named ``[PROCESSING]`` (see `example config file <https://github.com/SAR-ARD/ERS_NRB/blob/main/config.ini>`_).

Users might create several sections in the same configuration file with parameter values that correspond to different
processing scenarios (e.g., for different areas of interest). Note that each section must contain all necessary
configuration parameters even if only a few are varied between the sections.

Command Line Interface
----------------------
Once a configuration file has been created and all of its parameters have been properly defined, it can be used to start
the processor using the command line interface (CLI) tool provided with the ERS_NRB package.

The following options are currently available:

::

    ers_nrb --help

Print a help message for the CLI tool.

::

    ers_nrb -c /path/to/config.ini

Start the processor using parameters defined in the default section of a ``config.ini`` file.

::

    ers_nrb -c /path/to/config.ini -s SECTION_NAME

Start the processor using parameters defined in section ``SECTION_NAME`` of a ``config.ini`` file.
