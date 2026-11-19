from pathlib import Path
import json,warnings
import geopandas as g,numpy as np,shapely
from shapely.strtree import STRtree
from shapely.geometry import LineString
from shapely.ops import nearest_points
from convert_gcj02 import convert
out=Path('outputs/poi_anchored_v16');p=g.read_file(out/'poi_audit.gpkg');d=g.read_file(out/'blocks_projected.gpkg');roads=g.read_file('data/raw/osm/gis_osm_roads_free_1.shp').to_crs(32651);roads=roads[roads.fclass.isin(['primary','secondary','tertiary','residential','unclassified','pedestrian','footway','path','service'])&(roads.tunnel!='T')].copy();rg=roads.geometry.to_numpy();rt=STRtree(rg);bg=d.geometry.to_numpy();bt=STRtree(bg);rows=[];bd=[];rd=[]
for row in p.itertuples():
 point=row.geometry;near,dist=bt.query_nearest(point,max_distance=200,return_distance=True,all_matches=False);bd.append(float(dist[0]) if len(dist) else -1)
 near2,dist2=rt.query_nearest(point,max_distance=150,return_distance=True,all_matches=False);rd.append(float(dist2[0]) if len(dist2) else -1)
 if len(near2):
  idx=int(near2[0]);end=nearest_points(point,rg[idx])[1]
  if point.distance(end)>0.01:rows.append({'poi_id':row.poi_id,'kind':row.kind,'road_name':roads.name.iloc[idx] or '', 'road_class':roads.fclass.iloc[idx], 'distance_m':float(dist2[0]),'status':'nearest_road_not_verified_walk_route','geometry':LineString([point,end])})
p['land_dist_m']=bd;p['road_dist_m']=rd;p.to_file(out/'poi_audit_extended.gpkg',driver='GPKG')
for df,name in [(p,'poi_audit_extended_gcj02.geojson'),(g.GeoDataFrame(rows,crs=32651),'nearest_road_review_gcj02.geojson')]:
 if df.empty:continue
 x=df.to_crs(4326);x.geometry=shapely.transform(x.geometry.array,convert,include_z=False);x.crs=None
 with warnings.catch_warnings():
  warnings.filterwarnings('ignore',message="'crs' was not provided.*");x.to_file(out/name,driver='GeoJSON')
print(json.dumps({'poi':len(p),'land_within_200m':sum(v>=0 for v in bd),'road_within_150m':sum(v>=0 for v in rd),'links':len(rows)},indent=2))
