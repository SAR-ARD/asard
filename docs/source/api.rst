API Documentation
=================


Configuration
-------------

.. automodule:: ERS_NRB.config
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

.. automodule:: ERS_NRB.processor
    :members:
    :undoc-members:
    :show-inheritance:

SAR
^^^
`ERS_NRB` offers a mechanism to plug in different SAR processors.
The software offers a `snap` reference implementation module that can be translated to other processors.
All "main interface" functions need to be implemented so that `ERS_NRB` can fully interact with the module.

SNAP
++++

.. automodule:: ERS_NRB.snap
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

.. automodule:: ERS_NRB.ard
    :members:
    :undoc-members:
    :show-inheritance:

    .. autosummary::
        :nosignatures:

        append_metadata
        product_info

Ancillary Functions
-------------------

.. automodule:: ERS_NRB.ancillary
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

.. automodule:: ERS_NRB.metadata.extract
    :members:
    :undoc-members:
    :show-inheritance:

    .. autosummary::
        :nosignatures:

        get_prod_meta
        meta_dict
