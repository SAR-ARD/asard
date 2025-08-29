import pytest
import ERS_NRB.processor as process


@pytest.mark.parametrize(
    'sensor, acquisition_mode, date, aoi_tiles',
    [
        ('ASAR', 'APP', '2011-02-17', '31UFT'),
        ('ASAR', 'APS', '2003-07-15', '31UFT'),
        ('ASAR', 'IMP', '2011-08-20', '32UPB'),
        ('ASAR', 'IMS', '2004-07-03', '33TTG'),
    ]
)
def test_process(sensor, acquisition_mode, date, aoi_tiles,
                 testdata_dir, tmpdir, tmp_home):
    process.main(config_file=None,  # read processor default file
                 work_dir=str(tmpdir),
                 scene_dir=str(testdata_dir),
                 sensor=sensor,
                 mindate=date,
                 maxdate=date,
                 acq_mode=acquisition_mode,
                 mode='sar, nrb',
                 aoi_tiles=aoi_tiles,
                 gpt_args="-J-Xmx32G -c 22G -q 16"
                 )
