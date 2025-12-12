import pytest
import asard.processor as process


@pytest.mark.parametrize(
    'sensor, acquisition_mode, date, aoi_tiles',
    [
        ('ASAR', 'APP', '2002-12-21', '32TPS'),
        ('ASAR', 'APS', '2003-06-14', '32TPS'),
        ('ASAR', 'IMP', '2011-07-23', '32TPS'),
        ('ASAR', 'IMS', '2008-09-20', '32TPS'),
        ('ASAR', 'WSM', '2006-12-11', '32TPS'),
        ('ERS1', 'IMP', '1996-01-10', '32TPS'),
        ('ERS1', 'IMS', '1996-01-10', '32TPS'),
        ('ERS2', 'IMP', '1998-01-15', '32TPS'),
        ('ERS2', 'IMS', '1997-01-11', '32TPS')
    ]
)
@pytest.mark.parametrize("processor", ['snap', 'sarsenic'])
def test_process(processor, sensor, acquisition_mode, date, aoi_tiles,
                 testdata_dir, tmp_path, tmp_home):
    process.main(
        config_file=None,  # read processor default file
        processor=processor,
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
