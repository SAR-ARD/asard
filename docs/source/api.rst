API Documentation
=================


Configuration
-------------

.. automodule:: asard.config
    :members:
    :undoc-members:
    :show-inheritance:

    .. autosummary::
        :nosignatures:

        gdal_conf
        get_config
        get_keys
        read_config_file

Processing
----------

.. automodule:: asard.processor
    :members:
    :undoc-members:
    :show-inheritance:

SAR
^^^
`asard` offers a mechanism to plug in different SAR processors.
The software offers a `snap` reference implementation module that can be translated to other processors.
All "main interface" functions need to be implemented so that `asard` can fully interact with the module.

SNAP
++++

.. automodule:: asard.snap
    :members:
    :undoc-members:
    :show-inheritance:

    .. rubric:: main interface

    .. autosummary::
        :nosignatures:

        config_to_string
        get_config_keys
        get_config_section
        get_metadata


    .. rubric:: processor-specific functions

    .. autosummary::
        :nosignatures:

        process
        pre

ARD
^^^

.. automodule:: asard.ard
    :members:
    :undoc-members:
    :show-inheritance:

    .. autosummary::
        :nosignatures:

        append_metadata
        get_datasets
        format
        product_info

OSV file handling
-----------------

.. automodule:: asard.osv
    :members:
    :undoc-members:
    :show-inheritance:

    .. autosummary::
        :nosignatures:

        asar_osv_meta
        download_asar
        download_ers
        get

Ancillary Functions
-------------------

.. automodule:: asard.ancillary
    :members:
    :undoc-members:
    :show-inheritance:

    .. autosummary::
        :nosignatures:

        set_logging

Metadata
--------

Extraction
^^^^^^^^^^

.. automodule:: asard.metadata.extract
    :members:
    :undoc-members:
    :show-inheritance:

    .. autosummary::
        :nosignatures:

        get_prod_meta
        meta_dict
