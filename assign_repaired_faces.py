from pathlib import Path
import json
import numpy as np,geopandas as g,pandas as pd,shapely
from shapely.strtree import STRtree
out=Path('outputs/road_network_v2');blocks=g.read_file(out/'repaired_road_blocks.gpkg').to_crs(32651);old=g.read_file('outputs/revised_poc/activity_units_parent.gpkg').to_crs(32651);city=g.read_file('data/raw/shanghai_boundary.geojson').to_crs(32651).geometry.union_all();connectors=g.read_file(out/'inferred_gap_connectors.gpkg').to_crs(32651)
road_tree=STRtree(connectors.geometry.to_numpy());idx=STRtree(old.geometry.to_numpy());oldgeom=old.geometry.to_numpy();area=blocks.area.to_numpy();large=(~blocks.is_water.to_numpy())&(area>2e6);labels=[];counts={'water':0,'unresolved_large':0,'assigned':0,'unresolved_other':0};methods=[]
for i,r in enumerate(blocks.itertuples()):
 if r.is_water:labels.append('water');methods.append('water');counts['water']+=1;continue
 if large[i]:labels.append(f'unresolved-{i}');methods.append('unresolved_large');counts['unresolved_large']+=1;continue
 matches=idx.query(r.geometry,predicate='intersects')
 if not len(matches):labels.append(f'unresolved-{i}');methods.append('unresolved_other');counts['unresolved_other']+=1;continue
 # Assign a road-defined block to the previous candidate that covers the greatest area.
 overlap=shapely.area(shapely.intersection(oldgeom[matches],r.geometry));j=int(matches[np.argmax(overlap)])
 if overlap.max()<.01:labels.append(f'unresolved-{i}');methods.append('unresolved_other');counts['unresolved_other']+=1;continue
 labels.append(str(int(old.unit_seed.iloc[j])));methods.append('legacy_max_overlap');counts['assigned']+=1
blocks['assignment_key']=labels;blocks['assignment_method']=methods
print('face assignment',counts,flush=True)
# Any boundary touched by an inferred gap closure is flagged as such.
blocks['inferred_connector']=False
for i in range(len(blocks)):
 if blocks.is_water.iloc[i]:continue
 blocks.at[i,'inferred_connector']=len(road_tree.query(blocks.geometry.iloc[i].boundary,predicate='intersects'))>0
blocks.to_file(out/'road_based_blocks.gpkg',driver='GPKG')
land=blocks[~blocks.is_water].copy();groups=[];oldlookup=old.set_index('unit_seed')
for key,part in land.groupby('assignment_key',sort=False):
 geom=part.geometry.union_all()
 if key.startswith('unresolved-'):
  zone='road_network_incomplete';unit_seed=-int(key.split('-')[1])-1;inferred_crossing=False;anchor_source='none'
 else:
  unit_seed=int(key);v=oldlookup.loc[unit_seed];zone=v.zone_type;inferred_crossing=bool(v.inferred_crossing);anchor_source=v.anchor_source
 groups.append({'unit_key':key,'unit_seed':unit_seed,'zone_type':zone,'assignment_method':'legacy_max_overlap' if unit_seed>0 else 'unresolved_road_network','anchor_source':anchor_source,'inferred_crossing':inferred_crossing,'inferred_connector':bool(part.inferred_connector.any()),'geometry':geom})
units=g.GeoDataFrame(groups,crs=32651);units['area_km2']=units.area/1e6;water=blocks[blocks.is_water].copy();water['unit_key']='water';water['unit_seed']=0;water['zone_type']='water';water['assignment_method']='water';water=water[['unit_key','unit_seed','zone_type','assignment_method','geometry']]
partition=g.GeoDataFrame(pd.concat([units,water],ignore_index=True),crs=32651);print('units',len(units),'polygonizing outputs',flush=True)
units.to_file(out/'road_based_units.gpkg',driver='GPKG');partition.to_file(out/'road_based_partition.gpkg',driver='GPKG')
u=partition.geometry.union_all();info={'road_faces':len(blocks),'face_counts':counts,'land_units':len(units),'unresolved_large_faces':int(large.sum()),'unresolved_large_km2':float(area[large].sum()/1e6),'land_units_with_inferred_connectors':int(units.inferred_connector.sum()),'units_with_holes':int(sum(any(len(p.interiors) for p in (z.geoms if z.geom_type=='MultiPolygon' else [z])) for z in units.geometry)),'units_multipart':int((units.geom_type=='MultiPolygon').sum()),'invalid':int((~partition.is_valid).sum()),'gap_m2':float(city.difference(u).area),'overlap_m2':float(partition.area.sum()-u.area)};(out/'assignment_metrics.json').write_text(json.dumps(info,ensure_ascii=False,indent=2));print(json.dumps(info,ensure_ascii=False,indent=2),flush=True)
assert info['invalid']==0 and abs(info['overlap_m2'])<.1 and info['gap_m2']<1
