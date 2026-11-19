from pathlib import Path
import geopandas as g,pandas as pd,shapely,json
out=Path('outputs/revised_poc');out.mkdir(exist_ok=True)
city=g.read_file('data/raw/shanghai_boundary.geojson').to_crs(32651).geometry.union_all()
a=g.read_file('data/raw/osm/gis_osm_pois_a_free_1.shp').to_crs(32651);a=a[(a.fclass=='mall')&a.geometry.representative_point().within(city)].copy();a['geometry']=shapely.make_valid(a.geometry.array);a['source_type']='osm_mall_aoi';a['source_id']=a.osm_id.astype(str);a['center_wgs84_assumption']='WGS84';a=a[['source_id','name','source_type','geometry']]
c=g.read_file('outputs/city_imagery/city_commercial_core.geojson').to_crs(32651);c['source_id']=c.id.astype(str);c['name']='';c['source_type']='poi_density_core';c=c[['source_id','name','source_type','geometry']]
r=g.GeoDataFrame(pd.concat([a,c],ignore_index=True),crs=32651);r.to_crs(4326).to_file(out/'commercial_core_candidates.geojson',driver='GeoJSON');print('retained aoi',len(a),'density core',len(c),'total',len(r))
