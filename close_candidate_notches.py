from pathlib import Path
import json,warnings
import geopandas as g,numpy as np,shapely
from shapely.strtree import STRtree
from convert_gcj02 import convert
out=Path('outputs/notches_closed_v10');out.mkdir(exist_ok=True)
d=g.read_file('outputs/notch_example_v9/blocks_projected.gpkg').to_crs(32651);c=g.read_file('outputs/notch_example_v9/similar_notch_candidates.gpkg').to_crs(32651);geoms=d.geometry.to_numpy().copy();tree=STRtree(geoms);idmap={int(v):i for i,v in enumerate(d.piece_id)};water=g.read_file('data/raw/osm/gis_osm_water_a_free_1.shp').to_crs(32651);watergeom=water.geometry.to_numpy();wt=STRtree(watergeom);added=[];skipped=[]
# Candidates may overlap each other, so each is checked against prior additions.
for r in c.itertuples():
 i=idmap.get(int(r.piece_id));z=r.geometry
 if i is None:skipped.append((int(r.piece_id),'missing'));continue
 neighbors=[j for j in tree.query(z,predicate='intersects') if j!=i]
 if any(z.intersection(geoms[j]).area>.01 for j in neighbors):skipped.append((int(r.piece_id),'neighbor_overlap'));continue
 if any(z.intersection(watergeom[j]).area>.01 for j in wt.query(z,predicate='intersects')):skipped.append((int(r.piece_id),'water_overlap'));continue
 prior=[v['geometry'] for v in added if v['piece_id']!=int(r.piece_id) and v['geometry'].intersects(z)]
 if any(z.intersection(p).area>.01 for p in prior):skipped.append((int(r.piece_id),'candidate_overlap'));continue
 new=geoms[i].union(z)
 if new.geom_type=='MultiPolygon' and geoms[i].distance(z)<.001:
  from shapely.ops import nearest_points
  contact=nearest_points(geoms[i],z)[0];bridge=contact.buffer(.02);new=new.union(bridge)
  if any(bridge.intersection(geoms[j]).area>.01 for j in neighbors) or any(bridge.intersection(watergeom[j]).area>.01 for j in wt.query(bridge,predicate='intersects')):new=geoms[i].union(z)
 if new.geom_type!='Polygon' or not new.is_valid:skipped.append((int(r.piece_id),'geometry'));continue
 actual=new.difference(geoms[i]);geoms[i]=new
 if actual.area>.001:added.append({'piece_id':int(r.piece_id),'area_m2':actual.area,'geometry':actual})
d.geometry=geoms;d['area_km2']=d.area/1e6;d['holes_remaining']=[len(z.interiors) for z in geoms];d['notch_adjusted']=d.piece_id.isin([v['piece_id'] for v in added])|d.notch_adjusted.fillna(False)
d.to_file(out/'blocks_projected.gpkg',driver='GPKG');rev=g.GeoDataFrame(added,crs=32651);rev.to_file(out/'added_land.gpkg',driver='GPKG')
for frame,name in [(d,'blocks_gcj02.geojson'),(rev,'added_land_gcj02.geojson')]:
 x=frame.to_crs(4326);x.geometry=shapely.transform(x.geometry.array,convert,include_z=False);bad=~x.is_valid;x.loc[bad,'geometry']=shapely.make_valid(x.loc[bad,'geometry'].array);x=x.explode(index_parts=False,ignore_index=True);x=x[x.geom_type=='Polygon'];x.crs=None
 with warnings.catch_warnings():
  warnings.filterwarnings('ignore',message="'crs' was not provided.*");x.to_file(out/name,driver='GeoJSON')
m={'candidates':len(c),'closed':len(added),'skipped':len(skipped),'skip_reasons':{k:sum(s[1]==k for s in skipped) for k in set(s[1] for s in skipped)},'added_area_m2':sum(v['area_m2'] for v in added),'invalid':int((~d.is_valid).sum()),'skipped_piece_ids':skipped};(out/'metrics.json').write_text(json.dumps(m,ensure_ascii=False,indent=2));print(json.dumps({k:v for k,v in m.items() if k!='skipped_piece_ids'},indent=2))
