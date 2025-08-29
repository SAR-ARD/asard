Installation
============

SNAP
----

ERS_NRB requires ESA’s Sentinels Application Platform (SNAP) software to produce ERS_NRB products. It has been developed based on SNAP version 8.
Downloaders for different operating systems can be obtained from the `official webpage <https://step.esa.int/main/download/snap-download/>`_.

The following code can be used to replicate the software installation on a Linux OS:

::

    VERSION=8
    TARGET=~/SNAP"$VERSION"

    INSTALLER=esa-snap_sentinel_unix_"$VERSION"_0.sh
    wget https://download.esa.int/step/snap/"$VERSION".0/installers/"$INSTALLER"
    bash $INSTALLER -q -dir $TARGET
    $TARGET/bin/snap --nosplash --nogui --modules --update-all

See also the web page on how to `update SNAP from the command line <https://senbox.atlassian.net/wiki/spaces/SNAP/pages/30539785/Update+SNAP+from+the+command+line>`_.

Alternatively, updates for individual modules and versions can be downloaded in the `SNAP Update Center <https://step.esa.int/updatecenter/>`_.
The latest bundle that was used during release of version 0.1.2 is available here: https://step.esa.int/updatecenter/8.0_20220323-143356/.

ERS_NRB
-------

The ERS_NRB package is not yet available via conda-forge or other common package distribution channels. In the meantime,
the following shall provide a convenient installation option provided that Anaconda or Miniconda has been installed:

1. Create and then activate the conda environment

::

    conda env create --file https://raw.githubusercontent.com/SAR-ARD/ERS_NRB/main/environment-tpz.yaml
    conda activate nrb_env_tpz

2. Install the ERS_NRB package into the environment

The package version can be changed as necessary. See the `Tags <https://github.com/SAR-ARD/ERS_NRB/tags>`_ section of the
repository for available versions.

::

    pip install git+https://github.com/SAR-ARD/ERS_NRB.git@0.1.9

Docker
------

Both SNAP (9.0.0 + updates) and ERS_NRB can also be installed into a docker container using the Dockerfile that is provided with the package or using the dockerhub image in https://hub.docker.com/repository/docker/gr4n0t4/ers_nrb/general

