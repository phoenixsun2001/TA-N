from pathlib import Path
import json, os
import numpy as np, geopandas as gpd, shapely
from shapely import STRtree,make_valid
from shapely.ops import unary_union
import pandas as pd

root=Path('outputs/city_imagery');out=root/os.getenv('ROAD_FIT_OUTPUT','road_fit');orig=gpd.read_file(root/'closed/trade_areas_no_holes.geojson').to_crs(32651).sort_values('id').reset_index(drop=True);trial=gpd.read_file(out/'trial_candidates.gpkg').to_crs(32651).sort_values('id').reset_index(drop=True);zones=gpd.read_file(root/'closed/city_partition_no_holes.geojson').to_crs(32651);special=zones[zones.id<0].copy();city=gpd.read_file('data/raw/shanghai_boundary.geojson').to_crs(32651).geometry.union_all()
roads=gpd.read_file('data/raw/osm/gis_osm_roads_free_1.shp').to_crs(32651);roads=roads[roads.fclass.isin(['motorway','trunk','primary','secondary'])&(roads.bridge!='T')&(roads.tunnel!='T')&(roads.layer.fillna(0)==0)];roads=roads[roads.intersects(orig.geometry.union_all().buffer(500))].explode(index_parts=False);roads=roads[roads.geom_type=='LineString'];tree=STRtree(roads.geometry.to_numpy())
water=special[special.zone_type.isin(['osm_water_proxy','imagery_water_candidate'])].geometry.union_all()
centers=gpd.read_file('outputs/centers.geojson').to_crs(32651).set_index('id')
def close_length(geom):
 line=geom.boundary; length=line.length
 n=max(1,int(np.ceil(length/25))); points=shapely.line_interpolate_point(line,(np.arange(n)+.5)*length/n)
 idx=tree.query(points,predicate='dwithin',distance=30)
 return len(np.unique(idx[0]))*length/n
before=np.array([close_length(g) for g in orig.geometry]);after=np.array([close_length(g) for g in trial.geometry]);gain=after-before
keep=np.zeros(len(orig),bool);reason=[]
for i in range(len(orig)):
 a,b=orig.geometry[i],trial.geometry[i]; why=[]
 if not b.is_valid or b.geom_type!='Polygon' or len(b.interiors):why.append('invalid')
 if gain[i]<=0:why.append('no_improvement')
 if abs(b.area/a.area-1)>.10:why.append('area_change')
 if b.is_valid and b.difference(city).area>1:why.append('outside_city')
 if b.is_valid and b.difference(a).intersection(water).area>1:why.append('water_intrusion')
 if b.is_valid and not b.covers(centers.geometry.loc[orig.id[i]]):why.append('center_outside')
 keep[i]=not why;reason.append(','.join(why))
print('eligible',int(keep.sum()),flush=True)
# Revert lower-gain adjustments until all commercial interiors are disjoint.
while True:
 geom=np.array([trial.geometry[i] if keep[i] else orig.geometry[i] for i in range(len(orig))],dtype=object);idx=STRtree(geom);pairs=idx.query(geom,predicate='intersects');conf=[]
 for i,j in pairs.T:
  if i<j and geom[i].intersection(geom[j]).area>.01:conf.append((int(i),int(j)))
 if not conf:break
 i,j=conf[0];options=[x for x in (i,j) if keep[x]]
 if not options:raise RuntimeError('Original areas overlap; cannot safely fit')
 loser=min(options,key=lambda k:gain[k]);keep[loser]=False;reason[loser]='overlap_reverted'
print('accepted',int(keep.sum()),flush=True)
fit=orig.copy();fit.geometry=[trial.geometry[i] if keep[i] else orig.geometry[i] for i in range(len(orig))];fit['road_fit']=np.where(keep,'adjusted','unchanged');fit['close_before_m']=np.round(before,2);fit['close_after_m']=np.round(np.where(keep,after,before),2);fit['area_km2']=fit.area/1e6
fit_union=fit.geometry.union_all();original_union=orig.geometry.union_all();removed=original_union.difference(fit_union)
for idx,row in special.iterrows():
 geom=row.geometry
 if row.id==-1:geom=unary_union([geom,removed])
 special.at[idx,'geometry']=make_valid(geom.difference(fit_union))
special=special[~special.geometry.is_empty].copy();special['area_km2']=special.area/1e6
partition=gpd.GeoDataFrame(pd.concat([fit,special],ignore_index=True),crs=32651)
# WGS84 may turn microscopic artifacts into lines; keep polygonal components.
def polygon_only(geom):
 def parts(g):
  if g.geom_type=='Polygon':return [g]
  if hasattr(g,'geoms'):return [p for x in g.geoms for p in parts(x)]
  return []
 return unary_union(parts(make_valid(geom)))
wgs=partition.to_crs(4326);wgs.geometry=wgs.geometry.map(polygon_only);fw=fit.to_crs(4326)
assert all(fw.is_valid) and all(wgs.is_valid) and (fw.geom_type=='Polygon').all() and sum(len(g.interiors) for g in fw.geometry)==0
fw.to_file(out/'trade_areas_road_fit.geojson',driver='GeoJSON');wgs.to_file(out/'city_partition_road_fit.geojson',driver='GeoJSON')
roads.to_crs(4326).to_file(out/'candidate_major_roads.geojson',driver='GeoJSON')
# Verify exported coordinates in projected CRS.
check=gpd.read_file(out/'city_partition_road_fit.geojson').to_crs(32651);check_com=check[check.id>0];union=check.geometry.union_all();cross=check_com.area.sum()-check_com.geometry.union_all().area
metrics={'max_snap_m':float(os.getenv('ROAD_FIT_DISTANCE','150')),'sample_spacing_m':25,'closeness_m':30,'parallel_max_degrees':float(os.getenv('ROAD_FIT_ANGLE','40')),'road_types':['motorway','trunk','primary','secondary'],'selected_road_segments':len(roads),'trade_area_count':len(fit),'adjusted_count':int(keep.sum()),'unchanged_count':int((~keep).sum()),'eligible_count':int(sum(not x or x=='overlap_reverted' for x in reason)),'candidate_snapped_points':int(trial.snapped_points.sum()),'boundary_near_roads_before_pct':float(before.sum()/orig.length.sum()*100),'boundary_near_roads_after_pct':float(np.where(keep,after,before).sum()/fit.length.sum()*100),'mean_abs_area_change_pct_adjusted':float(np.mean(np.abs(fit.area[keep].to_numpy()/orig.area[keep].to_numpy()-1))*100) if keep.any() else 0,'commercial_overlap_m2':float(cross),'partition_overlap_m2':float(check.area.sum()-union.area),'partition_gap_m2':float(city.difference(union).area),'invalid_exported_count':int((~check.is_valid).sum()),'remaining_holes':int(sum(len(g.interiors) for g in check_com.geometry)),'rejection_reasons':{key:sum(key in x.split(',') for x in reason) for key in ['invalid','no_improvement','area_change','outside_city','water_intrusion','center_outside','overlap_reverted']}}
assert abs(metrics['commercial_overlap_m2'])<.1 and abs(metrics['partition_overlap_m2'])<.1 and metrics['partition_gap_m2']<.1 and metrics['invalid_exported_count']==0 and metrics['remaining_holes']==0
(out/'metrics.json').write_text(json.dumps(metrics,ensure_ascii=False,indent=2));print(json.dumps(metrics,ensure_ascii=False,indent=2))
