from pathlib import Path
import json
import geopandas as g,numpy as np,pandas as pd,shapely
from shapely.strtree import STRtree
out=Path('outputs/road_network_v2');b=g.read_file(out/'repaired_block_assignment.gpkg').to_crs(32651);seeds={x['seed_id']:x for x in json.loads((out/'road_seed_info.json').read_text())};connectors=g.read_file(out/'inferred_gap_connectors.gpkg').to_crs(32651);ct=STRtree(connectors.geometry.to_numpy());edges=g.read_file(out/'repaired_passage_links.gpkg');ids=b.unit_seed.to_numpy();risk=set()
for r in edges[edges.kind=='inferred_unverified'].itertuples():
 i,j=int(r.block_a),int(r.block_b)
 if ids[i]>0 and ids[i]==ids[j]:risk.add(int(ids[i]))
land=b[(~b.is_water)&(~b.large_flag)].copy();water=b[b.is_water].copy();large=b[(~b.is_water)&(b.large_flag)].copy();rows=[]
for sid,part in land.groupby('unit_seed'):
 if sid<=0:raise RuntimeError('Eligible land face lacks a seed')
 seed=seeds[int(sid)];labels=[x for x in seed['anchor_labels'] if x];geom=part.geometry.union_all()
 rows.append({'unit_seed':int(sid),'zone_type':'commercial_anchor' if seed['seed_kind'] not in ['local','local_component'] else 'local_activity','anchor_source':seed['seed_kind'],'anchor_count':len(labels),'anchor_names':' | '.join(labels)[:2000],'inferred_crossing':int(sid) in risk,'inferred_road_gap':bool(len(ct.query(geom.boundary,predicate='intersects'))),'geometry':geom})
for r in large.itertuples():
 rows.append({'unit_seed':-int(r.id)-1,'zone_type':'road_network_incomplete','anchor_source':'unresolved','anchor_count':0,'anchor_names':'','inferred_crossing':False,'inferred_road_gap':bool(len(ct.query(r.geometry.boundary,predicate='intersects'))),'geometry':r.geometry})
areas=g.GeoDataFrame(rows,crs=32651);areas['area_km2']=areas.area/1e6;water['unit_seed']=0;water['zone_type']='water';water=water[['unit_seed','zone_type','geometry']];part=g.GeoDataFrame(pd.concat([areas,water],ignore_index=True),crs=32651)
areas.to_file(out/'road_graph_units.gpkg',driver='GPKG');part.to_file(out/'road_graph_partition.gpkg',driver='GPKG')
city=g.read_file('data/raw/shanghai_boundary.geojson').to_crs(32651).geometry.union_all();u=part.geometry.union_all();m=json.loads((out/'road_graph_metrics.json').read_text());m.update({'land_units':len(areas),'commercial_units':int((areas.zone_type=='commercial_anchor').sum()),'local_units':int((areas.zone_type=='local_activity').sum()),'unresolved_units':int((areas.zone_type=='road_network_incomplete').sum()),'unresolved_km2':float(areas.area[areas.zone_type=='road_network_incomplete'].sum()/1e6),'units_with_inferred_road_gaps':int(areas.inferred_road_gap.sum()),'units_with_inferred_crossing':int(areas.inferred_crossing.sum()),'units_multipart':int((areas.geom_type=='MultiPolygon').sum()),'units_with_holes':int(sum(any(len(p.interiors) for p in (geom.geoms if geom.geom_type=='MultiPolygon' else [geom])) for geom in areas.geometry)),'invalid':int((~part.is_valid).sum()),'gap_m2':float(city.difference(u).area),'overlap_m2':float(part.area.sum()-u.area)})
(out/'road_graph_final_metrics.json').write_text(json.dumps(m,ensure_ascii=False,indent=2));print(json.dumps(m,ensure_ascii=False,indent=2),flush=True)
assert m['invalid']==0 and abs(m['overlap_m2'])<.1 and m['gap_m2']<1
