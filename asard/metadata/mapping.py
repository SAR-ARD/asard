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

# ERS-1 IM, ERS-2 IM and ENVISAT IM - Not Applicable
# ENVISAT AP and WS - Not Applied
NOISE_MAP = {'IMS': 'Not Applicable',
             'IMP': 'Not Applicable',
             'IMM': 'Not Applicable',
             'APP': 'Not Applied',
             'APS': 'Not Applied',
             'WSM': 'Not Applied',
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
                'title': 'Acquisition ID Image'},
    '-np-vv.tif': {'type': 'noise power VV',
                   'unit': None,
                   'role': 'noise-power',
                   'title': 'Noise Power VV'},
    '-np-vh.tif': {'type': 'noise power VH',
                   'unit': None,
                   'role': 'noise-power',
                   'title': 'Noise Power VH'},
    '-np-hh.tif': {'type': 'noise power HH',
                   'unit': None,
                   'role': 'noise-power',
                   'title': 'Noise Power HH'},
    '-np-hv.tif': {'type': 'noise power HV',
                   'unit': None,
                   'role': 'noise-power',
                   'title': 'Noise Power HV'}}

URL = {
    'ancillaryData_KML': 'https://sentiwiki.copernicus.eu/__attachments/1692737/'
                         'S2A_OPER_GIP_TILPAR_MPC__20151209T095117_V20150622T000000_21000101T000000_B00.zip',
    'card4l_nrb': 'https://ceos.org/ard/files/PFS/NRB/v5.5/CARD4L-PFS_NRB_v5.5.pdf',
    'card4l_orb': 'https://ceos.org/ard/files/PFS/ORB/v1.0/'
                  'CARD4L_Product_Family_Specification_Ocean_Radar_Backscatter-v1.0.pdf',
    'faradayRotationReference': None,
    'geoCorrAccuracyReference': None,
    'geoCorrAlgorithm': None,
    'griddingConventionURL': 'https://www.mgrs-data.org/data/documents/nga_mgrs_doc.pdf',
    'noiseRemovalAlgorithm': None,
    'orbitDataAccess': None,
    'radiometricAccuracyReference': None,
    'RTCAlgorithm': 'https://doi.org/10.1109/Tgrs.2011.2120616',
    'sensorCalibration': None,
    'source_access': 'https://dataspace.copernicus.eu',
    'source_doi': None,
    'windNormReferenceModel': None
}
