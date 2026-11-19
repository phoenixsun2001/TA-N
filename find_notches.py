from pathlib import Path
import geopandas as g,numpy as np,shapely,warnings
from shapely.strtree import STRtree
from convert_gcj02 import convert
out=Path('outputs/notch_example_v9');d=g.read_file('outputs/closed_small_holes_v8/blocks_projected_final.gpkg').to_crs(32651);a=d.geometry.to_numpy();tree=STRtree(a);w=g.read_file('data/raw/osm/gis_osm_water_a_free_1.shp').to_crs(32651);wg=w.geometry.to_numpy();wt=STRtree(wg);rows=[]
for i,p in enumerate(a):
 if p.area<20000:continue
 try: closed=p.buffer(8,quad_segs=4).buffer(-8,quad_segs=4).union(p);diff=closed.difference(p)
 except Exception:continue
 for part in shapely.get_parts(diff):
  if part.geom_type!='Polygon' or not 500<=part.area<=20000:continue
  bbox=part.bounds;long=max(bbox[2]-bbox[0],bbox[3]-bbox[1]);short=min(bbox[2]-bbox[0],bbox[3]-bbox[1]);
  if long<50 or long<short*2:continue
  near=tree.query(part,predicate='intersects');over=sum(part.intersection(a[j]).area for j in near if j!=i)
  if over>.1 or len(wt.query(part,predicate='intersects')):continue
  rows.append({'piece_id':int(d.piece_id.iloc[i]),'area_m2':part.area,'status':'candidate_unverified','geometry':part})
if rows:
 x=g.GeoDataFrame(rows,crs=32651);x.to_file(out/'similar_notch_candidates.gpkg',driver='GPKG');x=x.to_crs(4326);x.geometry=shapely.transform(x.geometry.array,convert,include_z=False);x.crs=None
 with warnings.catch_warnings():
  warnings.filterwarnings('ignore',message="'crs' was not provided.*");x.to_file(out/'similar_notch_candidates_gcj02.geojson',driver='GeoJSON')
print('candidates',len(rows),flush=True)
