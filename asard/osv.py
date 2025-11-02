import os
import re
import numpy as np
from urllib.parse import urlparse
from ftplib import FTP_TLS
from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED
import tarfile as tf
from datetime import datetime, timezone
from dateutil.relativedelta import relativedelta
from pyroSAR.drivers import ID
from pyroSAR.examine import ExamineSnap
from spatialist.ancillary import finder
import logging

log = logging.getLogger('asard')

osv_asar = (r'DOR_VOR_AXVF-P'
            r'(?P<publish>[0-9]{8}_[0-9]{6})_'
            r'(?P<start>[0-9]{8}_[0-9]{6})_'
            r'(?P<stop>[0-9]{8}_[0-9]{6})'
            r'(?:\.zip|)')


def download_asar(
        username: str,
        password: str,
        target: str | None = None,
        start: datetime | None = None,
        stop: datetime | None = None
) -> None:
    """
    Download DORIS OSV files from ftps://envorb-ftp-ds.eo.esa.int.
    
    
    Parameters
    ----------
    username:
        the FTP access username
    password:
        the FTP access password
    target:
        the download location. Default None: Store data in a subdirectory
        of the configured SNAP auxiliary data path
        (See :meth:`pyroSAR.examine.ExamineSnap.auxdatapath`).
    start:
        the acquisition start. Default None: no acquisition start limit
    stop:
        the acquisition stop. Default None: no acquisition stop limit

    Returns
    -------

    """
    if target is None:
        target = ExamineSnap().auxdatapath
    
    url = 'ftps://envorb-ftp-ds.eo.esa.int/Envisat_DORIS/DOR_VOR_AX_f'
    
    parsed = urlparse(url)
    
    ftp = FTP_TLS()
    ftp.connect(host=parsed.netloc, timeout=60, port=21)
    ftp.auth()
    ftp.login(username, password)
    ftp.prot_p()
    ftp.cwd(parsed.path)
    folders = ftp.nlst()
    
    if start is not None:
        if start.tzinfo is None:
            start = start.replace(tzinfo=timezone.utc)
        folder_min = start.strftime('%Y%m')
    else:
        folder_min = min(folders)
    
    if stop is not None:
        if stop.tzinfo is None:
            stop = stop.replace(tzinfo=timezone.utc)
        folder_max = stop.strftime('%Y%m')
    else:
        folder_max = max(folders)
    
    for folder in folders:
        if folder_min > folder or folder > folder_max:
            continue
        folder_path = parsed.path + '/' + folder
        ftp.cwd(folder_path)
        files = ftp.nlst()
        for file in files:
            meta = asar_osv_meta(file)
            if start is not None:
                if meta['stop'] < start:
                    continue
            if stop is not None:
                if meta['start'] > stop:
                    continue
            
            remote = f'{parsed.path}/{folder}/{file}'
            local = Path(target) / 'Orbits' / 'Doris' / 'vor' / folder / (file + '.zip')
            local.parent.mkdir(parents=True, exist_ok=True)
            if not local.exists():
                log.info(f'{local} <<-- {remote}')
                with ZipFile(local, mode="w", compression=ZIP_DEFLATED) as zf:
                    with zf.open(Path(remote).name, mode="w") as zf_file:
                        ftp.retrbinary(f"RETR {remote}", zf_file.write)
    ftp.quit()


def download_ers(
        username: str,
        password: str,
        satellite: str | None = None,
        target: str | None = None
) -> None:
    """
    Download REAPER OSV files from ftps://ersorb-ftp-ds.eo.esa.int.
    
    Parameters
    ----------
    username:
        the FTP access username
    password:
        the FTP access password
    satellite:
        Either ERS1, ERS2 or None (both satellites).
    target:
        the download location. Default None: Store data in a subdirectory
        of the configured SNAP auxiliary data path
        (See :meth:`pyroSAR.examine.ExamineSnap.auxdatapath`).

    Returns
    -------

    """
    satellite_options = ['ERS1', 'ERS2']
    if satellite is None:
        satellites = satellite_options
    else:
        if satellite not in satellite_options:
            raise ValueError(f'satellite {satellite} is not supported')
        else:
            satellites = [satellite]
    
    if target is None:
        target = ExamineSnap().auxdatapath
    
    url = 'ftps://ersorb-ftp-ds.eo.esa.int/ERS_Orbit_files/POD_REAPER_v2_FullMission'
    
    parsed = urlparse(url)
    
    ftp = FTP_TLS()
    ftp.connect(host=parsed.netloc, timeout=60, port=21)
    ftp.auth()
    ftp.login(username, password)
    ftp.prot_p()
    ftp.cwd(parsed.path)
    
    for sat in satellites:
        base = f'deos{sat.lower()}tar.gz'
        remote = f'{parsed.path}/{sat}/{base}'
        local = (Path(target) / 'Orbits' / 'REAPER' /
                 'POD_REAPER_v2_FullMission' / 'DEOS' / sat / base)
        local.parent.mkdir(parents=True, exist_ok=True)
        if not local.exists():
            log.info(f'{local} <<-- {remote}')
            with open(local, 'wb') as f:
                ftp.retrbinary(f"RETR {remote}", f.write)
            archive = tf.open(local, 'r')
            archive.extractall(local.parent)
            archive.close()
    ftp.quit()


def get(
        scene: ID,
        offline: bool = False,
        username: str | None = None,
        password: str | None = None,
        target: str | None = None
) -> str | None:
    """
    Get the matching OSV file for a satellite product.
    
    Parameters
    ----------
    scene:
        the satellite product object
    offline:
        only search offline?
    username:
        the FTP access username
    password:
        the FTP access password
    target:
        the download location. Default None: Store data in a subdirectory
        of the configured SNAP auxiliary data path
        (See :meth:`pyroSAR.examine.ExamineSnap.auxdatapath`).

    Returns
    -------
        either the name of an existing file or None if no file was found
    """
    if target is None:
        target = ExamineSnap().auxdatapath
    
    if scene.sensor == 'ASAR':
        if not offline:
            download_asar(username=username, password=password, target=target,
                          start=scene.start_dt, stop=scene.stop_dt)
        location = Path(target) / 'Orbits' / 'Doris' / 'vor'
        month_previous = (scene.start_dt - relativedelta(months=1)).strftime('%Y%m')
        month_current = scene.start_dt.strftime('%Y%m')
        folders = [location / month_previous, location / month_current]
        folders = [x for x in folders if x.exists()]
        candidates = []
        for folder in folders:
            files = finder(target=str(folder), matchlist=[osv_asar],
                           recursive=False, regex=True)
            for file in files:
                meta = asar_osv_meta(file)
                if meta['start'] < scene.start_dt < meta['stop']:
                    candidates.append((file, meta))
        if len(candidates) == 0:
            return None
        elif len(candidates) == 1:
            return candidates[0][0]
        else:
            diff_start = [scene.start_dt - x[1]['start'] for x in candidates]
            diff_stop = [x[1]['stop'] - scene.stop_dt for x in candidates]
            # determine the OSV file(s), whose range best covers that of the scene
            # a balanced score is computed optimizing the time range covered by the OSV
            # file before and after the scene acquisition
            a_num = np.array([x.total_seconds() for x in diff_start])
            b_num = np.array([x.total_seconds() for x in diff_stop])
            a_norm = (a_num - a_num.min() + 1) / (a_num.max() - a_num.min() + 1)
            b_norm = (b_num - b_num.min() + 1) / (b_num.max() - b_num.min() + 1)
            score = a_norm + b_norm
            indices = np.flatnonzero(score == score.max())
            if len(indices) == 1:
                return candidates[indices[0]][0]
            else:
                # if multiple files have the same score, return the last published file
                publish = [candidates[x][1]['publish'] for x in indices]
                idx = publish.index(max(publish))
                return candidates[idx][0]
    elif scene.sensor in ['ERS1', 'ERS2']:
        if not offline:
            download_ers(username=username, password=password,
                         target=target, satellite=scene.sensor)
        location = (Path(target) / 'Orbits' / 'REAPER' /
                    'POD_REAPER_v2_FullMission' / 'DEOS' / scene.sensor)
        prefix = 'xxo' if scene.sensor == 'ERS1' else 'adr'
        base = f'{prefix}.{scene.start[2:8]}.sp3'
        osv = location / base
        if not osv.exists():
            return None
        return str(osv)
    else:
        raise ValueError(f'unknown sensor: {scene.sensor}')


def asar_osv_meta(
        file: str
) -> dict[str, str | datetime]:
    """
    Read metadata from the ASAR OSV file name.
    Timestamps are converted to timezone-aware datetime objects.
    
    Parameters
    ----------
    file:
        the DORIS OSV file name.

    Returns
    -------
        the metadata dictionary
    """
    meta = re.match(osv_asar, os.path.basename(file)).groupdict()
    for date in ['publish', 'start', 'stop']:
        dt = datetime.strptime(meta[date], '%Y%m%d_%H%M%S')
        meta[date] = dt.replace(tzinfo=timezone.utc)
    return meta
