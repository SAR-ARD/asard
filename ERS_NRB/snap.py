import os
from spatialist.ancillary import finder
from pyroSAR.snap.auxil import parse_recipe


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
