from pathlib import Path
import json
import numpy as np,geopandas as g,pandas as pd,shapely
from shapely.ops import unary_union
out=Path('outputs/revised_poc');b=g.read_file(out/'activity_block_assignment.gpkg').to_crs(32651);active=b[b.active&~b.is_water];s=json.loads((out/'activity_seeds.json').read_text());lookup={v['seed_id']:v for v in s};links=g.read_file(out/'inferred_crossing_candidates.gpkg');ids=b.unit_seed.to_numpy();risk=set()
for r in links.itertuples():
 i,j=int(r.block_a),int(r.block_b)
 if ids[i]>0 and ids[i]==ids[j]:risk.add(int(ids[i]))
rows=[]
for k,part in active.groupby('unit_seed'):
 seed=lookup[int(k)];labels=[x for x in seed['anchor_labels'] if x]
 rows.append({'unit_seed':int(k),'zone_type':'commercial_anchor' if seed['seed_kind'] not in ['local','local_component'] else 'local_activity','anchor_source':seed['seed_kind'],'anchor_count':len(labels),'anchor_names':' | '.join(labels)[:2000],'fallback_split':bool(part.synthetic_split.any()),'inferred_crossing':int(k) in risk,'geometry':unary_union(part.geometry)})
areas=g.GeoDataFrame(rows,crs=32651);areas['area_km2']=areas.area/1e6;water=b[b.active&b.is_water].copy();water['unit_seed']=0;water['zone_type']='water';water=water[['unit_seed','zone_type','geometry']];part=g.GeoDataFrame(pd.concat([areas,water],ignore_index=True),crs=32651)
areas.to_file(out/'activity_units_parent.gpkg',driver='GPKG');part.to_file(out/'activity_partition.gpkg',driver='GPKG')
city=g.read_file('data/raw/shanghai_boundary.geojson').to_crs(32651).geometry.union_all();u=part.geometry.union_all();m=json.loads((out/'activity_metrics.json').read_text());m.update({'activity_units':len(areas),'water_components':len(water),'fallback_area_km2':float(active.area[active.synthetic_split].sum()/1e6),'land_area_km2':float(active.area.sum()/1e6),'units_with_fallback':int(areas.fallback_split.sum()),'units_with_inferred_crossing':int(areas.inferred_crossing.sum()),'invalid_geometries':int((~part.is_valid).sum()),'units_multipart':int((areas.geom_type=='MultiPolygon').sum()),'units_with_holes':sum(any(len(p.interiors) for p in (x.geoms if x.geom_type=='MultiPolygon' else [x])) for x in areas.geometry),'overlap_m2':float(part.area.sum()-u.area),'gap_m2':float(city.difference(u).area)})
(out/'activity_parent_metrics.json').write_text(json.dumps(m,ensure_ascii=False,indent=2));print(json.dumps(m,ensure_ascii=False,indent=2),flush=True)
assert m['invalid_geometries']==0 and abs(m['overlap_m2'])<.1 and m['gap_m2']<.1
