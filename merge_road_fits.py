from pathlib import Path
import json
import numpy as np, geopandas as gpd, pandas as pd, shapely
from shapely.ops import unary_union
from shapely import make_valid
root=Path('outputs/city_imagery');out=root/'road_fit';backup=root/'road_fit_75m'
base=gpd.read_file(root/'closed/trade_areas_no_holes.geojson').to_crs(32651).sort_values('id').reset_index(drop=True)
a=gpd.read_file(backup/'trade_areas_road_fit.geojson').to_crs(32651).sort_values('id').reset_index(drop=True)
b=gpd.read_file(out/'trade_areas_road_fit.geojson').to_crs(32651).sort_values('id').reset_index(drop=True)
selected=a.copy(); selected['fit_distance']=np.where(selected.road_fit=='adjusted',75,0)
# Try the 150 m candidates where the 75 m proposal was not accepted, without displacing any accepted neighbour.
for i in range(len(selected)):
 if selected.fit_distance[i] or b.road_fit[i]!='adjusted':continue
 geom=b.geometry[i];others=selected.drop(index=i);ids=others.sindex.query(geom,predicate='intersects')
 if any(geom.intersection(others.geometry.iloc[j]).area>.01 for j in ids):continue
 selected.at[i,'geometry']=geom;selected.at[i,'road_fit']='adjusted';selected.at[i,'fit_distance']=150;selected.at[i,'close_after_m']=b.close_after_m[i]
selected['area_km2']=selected.area/1e6
zones=gpd.read_file(root/'closed/city_partition_no_holes.geojson').to_crs(32651);special=zones[zones.id<0].copy();city=gpd.read_file('data/raw/shanghai_boundary.geojson').to_crs(32651).geometry.union_all();fit_union=selected.geometry.union_all();removed=base.geometry.union_all().difference(fit_union)
for idx,row in special.iterrows():
 geom=row.geometry
 if row.id==-1:geom=unary_union([geom,removed])
 special.at[idx,'geometry']=make_valid(geom.difference(fit_union))
special=special[~special.geometry.is_empty].copy();special['area_km2']=special.area/1e6
part=gpd.GeoDataFrame(pd.concat([selected,special],ignore_index=True),crs=32651)
def polygon_only(geom):
 def parts(g):
  if g.geom_type=='Polygon':return [g]
  if hasattr(g,'geoms'):return [p for x in g.geoms for p in parts(x)]
  return []
 return unary_union(parts(make_valid(geom)))
fw=selected.to_crs(4326);pw=part.to_crs(4326);pw.geometry=pw.geometry.map(polygon_only);assert fw.is_valid.all() and pw.is_valid.all()
fw.to_file(out/'trade_areas_road_fit.geojson',driver='GeoJSON');pw.to_file(out/'city_partition_road_fit.geojson',driver='GeoJSON')
check=gpd.read_file(out/'city_partition_road_fit.geojson').to_crs(32651);c=check[check.id>0];u=check.geometry.union_all();old=json.loads((out/'metrics.json').read_text());old.update({'method':'combined_75m_then_150m','adjusted_count':int((selected.fit_distance>0).sum()),'adjusted_75m':int((selected.fit_distance==75).sum()),'adjusted_150m':int((selected.fit_distance==150).sum()),'unchanged_count':int((selected.fit_distance==0).sum()),'boundary_near_roads_before_pct':float(selected.close_before_m.sum()/base.length.sum()*100),'boundary_near_roads_after_pct':float(selected.close_after_m.sum()/selected.length.sum()*100),'mean_abs_area_change_pct_adjusted':float(np.mean(np.abs(selected.area[selected.fit_distance>0].to_numpy()/base.area[selected.fit_distance>0].to_numpy()-1))*100),'commercial_overlap_m2':float(c.area.sum()-c.geometry.union_all().area),'partition_overlap_m2':float(check.area.sum()-u.area),'partition_gap_m2':float(city.difference(u).area),'invalid_exported_count':int((~check.is_valid).sum()),'remaining_holes':int(sum(len(g.interiors) for g in c.geometry))})
for k in ['eligible_count','candidate_snapped_points','rejection_reasons','parallel_max_degrees']: old.pop(k,None)
old['candidate_settings']=[{'max_snap_m':75,'parallel_max_degrees':30},{'max_snap_m':150,'parallel_max_degrees':40}]
assert abs(old['commercial_overlap_m2'])<.1 and abs(old['partition_overlap_m2'])<.1 and old['partition_gap_m2']<.1 and old['invalid_exported_count']==0 and old['remaining_holes']==0
(out/'metrics.json').write_text(json.dumps(old,ensure_ascii=False,indent=2));print(json.dumps(old,ensure_ascii=False,indent=2))
