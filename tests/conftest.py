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
        'ASA_APP_1PNESA20021221_205441_000000182012_00172_04234_0000.N1',
        'ASA_APS_1PNESA20030614_205444_000000172017_00172_06739_0000.N1',
        'ASA_IMP_1PRDAM20110723_211202_000000183104_00417_49138_0000.N1',
        'ASA_IMS_1PNESA20080920_205443_000000172072_00172_34294_0000.N1',
        'ASA_WSM_1PNDSI20061211_093150_000004532053_00394_24997_0000.N1',
        'SAR_IMP_1PRDAM19951206_100756_00000018G151_00437_22970_0000.E1',
        'SAR_IMS_1PNESA19960110_100754_00000018G152_00437_23471_0000.E1',
        'SAR_IMP_1PNESA19980115_100749_00000018A028_00437_14319_0000.E2',
        'SAR_IMS_1PNESA19970111_212446_00000018A018_00172_09044_0000.E2'
    ]
    
    for scene in scenes:
        fname = out / scene
        if not fname.is_file():
            raise RuntimeError(f'missing test data: {fname}')
    
    return out
