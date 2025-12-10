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
        out = Path(os.environ['ASARD_TESTDATA'])
    except KeyError:
        raise RuntimeError('ASARD_TESTDATA environment variable not set')
    if not out.is_dir():
        raise RuntimeError(f'ASARD_TESTDATA is not an existing directory: {out}')
    
    scenes = [
        'ASA_APP_1PNESA20110217_101056_000000163099_00324_46890_0000.N1',
        'ASA_APS_1PNESA20030715_095933_000000172018_00108_07176_0000.N1',
        'ASA_IMP_1PNESA20110820_092721_000000173105_00381_49533_0000.N1',
        'ASA_IMS_1PNESA20080909_140200_000000162072_00010_34132_0000.N1',
        'SAR_IMP_1PNESA19960301_100148_00000018G154_00165_24201_0000.E1',
        'SAR_IMS_1PNESA19940624_171232_00000017E140_00457_15378_0000.E1',
        'SAR_IMP_1PNESA19970401_060315_00000018A020_00306_10180_0000.E2',
        'SAR_IMS_1PNESA20021010_153746_00000016A078_00140_39072_0000.E2'
    ]
    
    for scene in scenes:
        fname = out / scene
        if not fname.is_file():
            raise RuntimeError(f'missing test data: {fname}')
    
    return out
