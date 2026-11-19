from pathlib import Path
import json,collections
import numpy as np,geopandas as g,pandas as pd,shapely
from shapely.geometry import Polygon,LineString
from shapely.ops import unary_union,split,triangulate
out=Path('outputs/revised_poc');b=g.read_file(out/'activity_block_assignment.gpkg').to_crs(32651);active=b[b.active&~b.is_water].copy();seeds=json.loads((out/'activity_seeds.json').read_text());seedmap={r['seed_id']:r for r in seeds};inferred=g.read_file(out/'inferred_crossing_candidates.gpkg');ids=b.unit_seed.to_numpy();risky=set()
for r in inferred.itertuples():
 i,j=int(r.block_a),int(r.block_b)
 if ids[i]>0 and ids[i]==ids[j]:risky.add(int(ids[i]))
# Split across an interior hole instead of filling it and swallowing other activity units.
def polygon_parts(geom):
 if geom.geom_type=='Polygon':return [geom]
 if hasattr(geom,'geoms'):return [p for v in geom.geoms for p in polygon_parts(v)]
 return []
def remove_holes(p):
 todo=[p];done=[];steps=0
 while todo:
  q=todo.pop();steps+=1
  if steps>10000:raise RuntimeError('hole split limit exceeded')
  if not q.interiors:done.append(q);continue
  hole=Polygon(q.interiors[0]);pt=hole.representative_point();minx,miny,maxx,maxy=q.bounds
  choices=[LineString([(minx-1,pt.y),(maxx+1,pt.y)]),LineString([(pt.x,miny-1),(pt.x,maxy+1)])];pieces=None
  for line in choices:
   candidates=[z for z in polygon_parts(split(q,line)) if z.area>1e-9]
   if len(candidates)>1 and sum(len(z.interiors) for z in candidates)<len(q.interiors):pieces=candidates;break
  if pieces is None:
   pieces=[]
   for t in triangulate(q):pieces.extend(z for z in polygon_parts(t.intersection(q)) if z.area>1e-9)
   if any(z.interiors for z in pieces):raise RuntimeError('hole-free decomposition failed')
  todo.extend(pieces)
 return done
records=[];holes=0;before_parts=0
for k,grp in active.groupby('unit_seed',sort=True):
 seed=seedmap[int(k)];merged=unary_union(grp.geometry);polys=polygon_parts(merged);before_parts+=len(polys);holes+=sum(len(p.interiors) for p in polys);pieces=[]
 for p in polys:pieces.extend(remove_holes(p))
 for j,p in enumerate(pieces,1):
  records.append({'activity_id':f'{int(k)}-{j:03d}','unit_seed':int(k),'unit_type':'commercial_anchor' if seed['seed_kind'] not in ['local','local_component'] else 'local_activity','anchor_source':seed['seed_kind'],'anchor_count':len([v for v in seed['anchor_labels'] if v]),'anchor_names':' | '.join([v for v in seed['anchor_labels'] if v])[:2000],'fallback_split':bool(grp.synthetic_split.any()),'inferred_crossing':int(k) in risky,'geometry':p})
print('units',len(seedmap),'polygons',len(records),'rings split',holes,flush=True)
areas=g.GeoDataFrame(records,crs=32651);areas['area_km2']=areas.area/1e6;water=b[b.active&b.is_water].geometry.union_all();city=g.read_file('data/raw/shanghai_boundary.geojson').to_crs(32651).geometry.union_all();land_union=areas.geometry.union_all();water=water.intersection(city).difference(land_union);water_gdf=g.GeoDataFrame([{'activity_id':'water','unit_seed':0,'unit_type':'water','geometry':water}],crs=32651)
part=g.GeoDataFrame(pd.concat([areas,water_gdf],ignore_index=True),crs=32651);part['zone_type']=part.unit_type
# Repair collapsed microscopic slivers after reprojection, retaining only polygonal components.
def clean(geom):return unary_union(polygon_parts(shapely.make_valid(geom)))
aw=areas.to_crs(4326);pw=part.to_crs(4326);invalid=~pw.is_valid;pw.loc[invalid,'geometry']=pw.loc[invalid,'geometry'].map(clean);assert aw.is_valid.all() and (aw.geom_type=='Polygon').all() and all(not len(x.interiors) for x in aw.geometry) and pw.is_valid.all()
aw.to_file(out/'activity_units.geojson',driver='GeoJSON');pw.to_file(out/'activity_partition.geojson',driver='GeoJSON')
check=g.read_file(out/'activity_partition.geojson').to_crs(32651);land=check[check.unit_seed>0];union=check.geometry.union_all();metrics=json.loads((out/'activity_metrics.json').read_text());metrics.update({'output_polygon_components':len(areas),'rings_cut':holes,'fallback_area_km2':float(active.area[active.synthetic_split].sum()/1e6),'land_area_km2':float(active.area.sum()/1e6),'areas_using_fallback':int(areas.fallback_split.sum()),'areas_with_inferred_crossing':int(areas.inferred_crossing.sum()),'invalid_geometries':int((~check.is_valid).sum()),'activity_holes':sum(len(x.interiors) for x in land.geometry),'activity_nonpolygon':int((land.geom_type!='Polygon').sum()),'overlap_m2':float(check.area.sum()-union.area),'gap_m2':float(city.difference(union).area)})
(out/'activity_final_metrics.json').write_text(json.dumps(metrics,ensure_ascii=False,indent=2));print(json.dumps(metrics,ensure_ascii=False,indent=2),flush=True)
assert metrics['invalid_geometries']==metrics['activity_holes']==metrics['activity_nonpolygon']==0 and abs(metrics['overlap_m2'])<.1 and metrics['gap_m2']<.1
