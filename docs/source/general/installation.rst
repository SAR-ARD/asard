Installation
============

asard
-----

The asard package is not yet available via conda-forge or other common package distribution channels. In the meantime,
the following shall provide a convenient installation option provided that Anaconda or Miniconda has been installed:

::

    mamba env create --file https://raw.githubusercontent.com/SAR-ARD/asard/main/environment.yaml
    mamba activate asard
    pip install --no-deps git+https://github.com/SAR-ARD/asard.git

SAR Processors
--------------

asard currently offers two different SAR processors: sarsenic and SNAP.

SARSENIC
^^^^^^^^

This processor uses an extension of sarsen to process Envisat and ERS images.
The name sarsenic was chosen because the asard module cannot have the same name as the imported sarsen software.
Two more packages need to be installed: asar-xarray for reading the products into xarray and envisat_sarsen, which is a fork of sarsen.
The previously created asard environment is expected to be activated.

::

    mamba env update --name asard --file https://raw.githubusercontent.com/SAR-ARD/asar-xarray/main/environment.yml
    pip install --no-deps git+https://github.com/SAR-ARD/asar-xarray.git

::

    mamba env update --name asard --file https://raw.githubusercontent.com/SAR-ARD/envisat_sarsen/main/environment.yml
    pip install --no-deps git+https://github.com/SAR-ARD/envisat_sarsen.git

SNAP
^^^^

The second processor is ESA’s Sentinels Application Platform (SNAP). asard has been developed based on SNAP version 12.
Downloaders for different operating systems can be obtained from the `official webpage <https://step.esa.int/main/download/snap-download/>`_.

The following code can be used to replicate the software installation on a Linux OS:

::

    VERSION=12
    TARGET=~/SNAP"$VERSION"

    INSTALLER=esa-snap_all_linux-"$VERSION".0.0.sh
    wget https://download.esa.int/step/snap/"$VERSION".0/installers/"$INSTALLER"
    bash $INSTALLER -q -dir $TARGET
    $TARGET/bin/snap --nosplash --nogui --modules --update-all

See also the web page on how to `update SNAP from the command line <https://senbox.atlassian.net/wiki/spaces/SNAP/pages/30539785/Update+SNAP+from+the+command+line>`_.

Alternatively, updates for individual modules and versions can be downloaded in the `SNAP Update Center <https://step.esa.int/updatecenter/>`_.
