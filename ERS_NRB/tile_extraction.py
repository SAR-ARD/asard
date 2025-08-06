from spatialist.vector import Vector
from shapely.geometry import box
from shapely.affinity import scale
from pyproj import CRS
from pyproj.aoi import AreaOfInterest
from pyproj.database import query_utm_crs_info
from pyproj import Transformer

from s1ard.tile_extraction import aoi_from_tile, tile_from_aoi


def main(config, tr):
    """
    
    Parameters
    ----------
    config: dict
        Dictionary of the parsed config parameters for the current process.
    tr: int
        The target pixel spacing in meters, which is passed to `pyroSAR.snap.util.geocode`.
    
    Returns
    -------
    geo_dict: dict
        Dictionary containing geospatial information for each unique MGRS tile ID that will be processed.
    align_dict: dict
        Dictionary containing 'xmax' and 'ymin' coordinates to use for the `standardGridOriginX` and
        `standardGridOriginY` parameters of `pyroSAR.snap.util.geocode`
    """
    
    if config['processing']['aoi_tiles'] is not None:
        tiles = config['processing']['aoi_tiles']
    elif config['processing']['aoi_tiles'] is None and config['processing']['aoi_geometry'] is not None:
        with Vector(config['processing']['aoi_geometry']) as aoi:
            tiles = tile_from_aoi(aoi)
    else:
        raise RuntimeError("Either 'aoi_tiles' or 'aoi_geometry' need to be provided!")
    
    geo_dict = {}
    for tile in tiles:
        with aoi_from_tile(tile=tile) as vec:
            ext = vec.extent
            epsg = vec.getProjection('epsg')
            xmax = ext['xmax'] - tr / 2
            ymin = ext['ymin'] + tr / 2
            geo_dict[tile] = {'ext': ext,
                              'epsg': epsg,
                              'xmax': xmax,
                              'ymin': ymin}
    
    align_dict = {'xmax': max([geo_dict[tile]['xmax'] for tile in list(geo_dict.keys())]),
                  'ymin': min([geo_dict[tile]['ymin'] for tile in list(geo_dict.keys())])}
    
    return geo_dict, align_dict


def no_aoi(ids, spacing):
    """
    
    Parameters
    ----------
    ids: list of pyroSAR.drivers.ID
        List containing the scenes to process.
    spacing: int
        The target pixel spacing in meters, which is passed to `pyroSAR.snap.util.geocode`.
    
    Returns
    -------
    geo_dict: dict
        Dictionary containing geospatial information for each scene that will be processed.
    align_dict: dict
        Dictionary containing 'xmax' and 'ymin' coordinates to use for the `standardGridOriginX` and
        `standardGridOriginY` parameters of `pyroSAR.snap.util.geocode`
    """
    
    geo_dict = {}
    for i, scene in enumerate(ids):
        # create transformer to convert points from lat lon to a UTM
        corners = scene.getCorners()
        utm_crs_list = query_utm_crs_info(
            datum_name='WGS 84',
            area_of_interest=AreaOfInterest(
                west_lon_degree=corners['xmin'],
                south_lat_degree=corners['ymin'],
                east_lon_degree=corners['xmin'],
                north_lat_degree=corners['ymin']))
        utm_crs = CRS.from_epsg(utm_crs_list[0].code)
        crs_lat_lon = CRS.from_epsg(4326)
        ## transform the LA tLON points to UTM
        transformer = Transformer.from_crs(crs_lat_lon, utm_crs, always_xy=True)
        ymin_utm, xmin_utm = transformer.transform(corners['xmin'], corners['ymin'])
        ymax_utm, xmax_utm = transformer.transform(corners['xmax'], corners['ymax'])
        
        ## create polygon using UTM points
        UTM_poly = box(ymin_utm, xmin_utm, ymax_utm, xmax_utm)
        UTM_poly_grown = scale(UTM_poly, xfact=1.2, yfact=1.2, zfact=0.0,
                               origin='center')  # grow the polygon by 20 % to ensure the AOI is covered
        
        ext = {'xmin': UTM_poly_grown.bounds[0], 'xmax': UTM_poly_grown.bounds[2], 'ymin': UTM_poly_grown.bounds[1],
               'ymax': UTM_poly_grown.bounds[3]}
        xmax = ext['xmax'] - spacing / 2
        ymin = ext['ymin'] + spacing / 2
        geo_dict[scene.outname_base()] = {'ext': ext,
                                          'epsg': int(utm_crs_list[0].code),
                                          'xmax': xmax,
                                          'ymin': ymin}
        align_dict = {'xmax': max([geo_dict[tile]['xmax'] for tile in list(geo_dict.keys())]),
                      'ymin': min([geo_dict[tile]['ymin'] for tile in list(geo_dict.keys())])}
    
    return geo_dict, align_dict
