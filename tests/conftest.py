import os
import pytest
import platform
from pathlib import Path


@pytest.fixture
def tmp_home(monkeypatch, tmp_path):
    home = tmp_path / 'tmp_home'
    home.mkdir()
    var = 'USERPROFILE' if platform.system() == 'Windows' else 'HOME'
    monkeypatch.setenv(var, str(home))
    assert os.path.expanduser('~') == str(home)
    yield home


@pytest.fixture
def testdata_dir():
    try:
        out = Path(os.environ['ERS_NRB_TESTDATA'])
    except KeyError:
        raise RuntimeError('ERS_NRB_TESTDATA environment variable not set')
    if not out.is_dir():
        raise RuntimeError(f'ERS_NRB_TESTDATA is not an existing directory: {out}')
    
    scenes = [
        'ASA_APP_1PNESA20110217_101056_000000163099_00324_46890_0000.N1',
        'ASA_APS_1PNESA20030715_095933_000000172018_00108_07176_0000.N1',
        'ASA_IMP_1PNESA20110820_092721_000000173105_00381_49533_0000.N1',
        'ASA_IMS_1PNESA20040703_205338_000000182028_00172_12250_00001672562030318361237.N1'
    ]
    
    for scene in scenes:
        fname = out / scene
        if not fname.is_file():
            raise RuntimeError(f'missing test data: {fname}')
    
    return out
