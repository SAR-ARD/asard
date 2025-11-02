Installation
============

SNAP
----

asard requires ESA’s Sentinels Application Platform (SNAP) software to produce ERS/ASAR NRB products. It has been developed based on SNAP version 12.
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

asard
-------

The asard package is not yet available via conda-forge or other common package distribution channels. In the meantime,
the following shall provide a convenient installation option provided that Anaconda or Miniconda has been installed:

1. Create and then activate the conda environment

::

    conda env create --file https://raw.githubusercontent.com/SAR-ARD/asard/main/environment.yaml
    conda activate nrb_env

2. Install the asard package into the environment

The package version can be changed as necessary. See the `Tags <https://github.com/SAR-ARD/asard/tags>`_ section of the
repository for available versions.

::

    pip install git+https://github.com/SAR-ARD/asard.git@0.1.9
