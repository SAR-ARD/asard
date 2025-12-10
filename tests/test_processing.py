import pytest
import asard.processor as process


@pytest.mark.parametrize(
    'sensor, acquisition_mode, date, aoi_tiles',
    [
        ('ASAR', 'APP', '2011-02-17', '31UFT'),
        ('ASAR', 'APS', '2003-07-15', '31UFT'),
        ('ASAR', 'IMP', '2011-08-20', '32UPB'),
        ('ASAR', 'IMS', '2008-09-09', '19JFN'),
        ('ASAR', 'WSM', '2006-12-11', '32TPS'),
        ('ERS1', 'IMP', '1996-03-01', '33VWE'),
        ('ERS1', 'IMS', '1994-06-24', '43RGL'),
        ('ERS2', 'IMP', '1997-04-01', '11TME'),
        ('ERS2', 'IMS', '2002-10-10', '47QQG')
    ]
)
@pytest.mark.parametrize("processor", ['snap', 'sarsenic'])
def test_process(processor, sensor, acquisition_mode, date, aoi_tiles,
                 testdata_dir, tmp_path, tmp_home):
    process.main(config_file=None,  # read processor default file
                 work_dir=str(tmp_path),
                 scene_dir=str(testdata_dir),
                 sensor=sensor,
                 mindate=date,
                 maxdate=date,
                 acq_mode=acquisition_mode,
                 mode='sar, nrb',
                 aoi_tiles=aoi_tiles,
                 gpt_args="-J-Xmx32G -c 22G -q 16"
                 )
