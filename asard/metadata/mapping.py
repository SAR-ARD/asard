from cesard.metadata.mapping import URL

ARD_PATTERN = r'^(?P<mission>ER[12]|ENV)' \
              r'(?P<sensor>A|S)' \
              r'(?P<mode>IM|AP|WS)' \
              r'(?P<product>NRB)_' \
              r'(?P<start>[0-9]{8}T[0-9]{6})_' \
              r'(?P<duration>[0-9]{4})_' \
              r'(?P<platform>_)' \
              r'(?P<rel_orbit>[0-9A-Z]{3})_' \
              r'(?P<aux_data_level>S)' \
              r'(?P<id>[0-9A-Z]{3})_' \
              r'(?P<phase>[0-9A-Z]{1})' \
              r'(?P<cycle>[0-9]{3})_' \
              r'(?P<pols>DH|DV|SH|SV)'

# Envisat
# FP = FOS predicted orbit state vectors (NRT processing)
# DN = DORIS Level 0 navigator product acquired at PDHS (NRT)
# FR = FOS restituted orbit state vectors
# DI = DORIS initial (preliminary) orbit
# DP = DORIS precise orbit If not used, set to ØØ.    
ORB_MAP = {'PD': 'predicted',
           'RS': 'restituted',
           'PC': 'precise',
           'PL': 'preliminary',
           'FP': 'predicted',
           'DN': 'navigator',
           'FR': 'restituted',
           'DI': 'preliminary',
           'DP': 'precise'}

# IM* - Not Applicable
# AP*, WS* - Not Applied
NOISE_MAP = {
    'IMS': False,
    'IMP': False,
    'IMM': False,
    'APP': False,
    'APS': False,
    'WSM': False,
    'WSS': False
}

SAMPLE_MAP = {
    '-dm.tif': {'type': 'mask',
                'unit': 'mask',
                'role': 'data-mask',
                'title': 'Data Mask Image',
                'values': {0: 'not layover, nor shadow',
                           1: 'layover',
                           2: 'shadow',
                           3: 'layover and shadow',
                           4: 'ocean water'}},
    '-ei.tif': {'type': 'angle',
                'unit': 'deg',
                'role': 'ellipsoid-incidence-angle',
                'title': 'Ellipsoid Incidence Angle'},
    '-lc.tif': {'type': 'scattering area',
                'unit': 'square_meters',
                'role': 'contributing-area',
                'title': 'Local Contributing Area'},
    '-li.tif': {'type': 'angle',
                'unit': 'deg',
                'role': 'local-incidence-angle',
                'title': 'Local Incidence Angle'},
    '-gs.tif': {'type': 'ratio',
                'unit': 'ratio',
                'role': 'gamma-sigma-ratio',
                'title': 'Gamma0 RTC to sigma0 RTC ratio'},
    '-id.tif': {'type': 'mask',
                'unit': None,
                'role': 'acquisition-id',
                'title': 'Acquisition ID Image'}
}

URL.update({
    'faradayRotationReference': None,
    'geoCorrAccuracyReference': None,
    'geoCorrAlgorithm': 'TBD',
    'noiseRemovalAlgorithm': 'TBD',
    'orbitDataAccess': 'TBD',
    'radiometricAccuracyReference': None,
    'RTCAlgorithm': 'https://doi.org/10.1109/Tgrs.2011.2120616',
    'sensorCalibration': None,
    'source_access': 'https://esar-ds.eo.esa.int',
    'source_doi': None
})
