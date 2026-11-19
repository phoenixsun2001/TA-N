import geopandas as g, numpy as np, shapely,json,warnings
from shapely.geometry import Polygon
p='outputs/riverfront_v15/blocks_gcj02.geojson';d=g.read_file(p);old=d.to_crs(32651);holes=sum(len(p.interiors) for p in d.geometry);d.geometry=[Polygon(p.exterior) for p in d.geometry];added=d.to_crs(32651).area.sum()-old.area.sum();d.crs=None
with warnings.catch_warnings():
 warnings.filterwarnings('ignore',message="'crs' was not provided.*");d.to_file('outputs/riverfront_v15/blocks_gcj02_exterior.geojson',driver='GeoJSON')
print(json.dumps({'holes_closed':holes,'pseudo_area_added_km2':added/1e6,'remaining_holes':sum(len(z.interiors) for z in d.geometry),'invalid':int((~d.is_valid).sum())},indent=2))
