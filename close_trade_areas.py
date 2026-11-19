from pathlib import Path
import json
import geopandas as gpd
from shapely.geometry import Polygon
from shapely import make_valid

src=Path('outputs/city_imagery');out=src/'closed';out.mkdir(exist_ok=True)
z=gpd.read_file(src/'city_fused.geojson').to_crs(32651)
commercial=z[z.id>0].copy();special=z[z.id<0].copy();boundary=gpd.read_file('data/raw/shanghai_boundary.geojson').to_crs(32651).geometry.union_all()
assert len(commercial)==190 and commercial.geometry.geom_type.eq('Polygon').all()
records=[]
for row in commercial.itertuples():
 for n,ring in enumerate(row.geometry.interiors,1):
  hole=Polygon(ring)
  intersect_types=[]
  for other in special.itertuples():
   area=hole.intersection(other.geometry).area
   if area>0.001:intersect_types.append({'zone_type':other.zone_type,'area_m2':area})
  records.append({'id':int(row.id),'hole_number':n,'area_m2':hole.area,'contents':json.dumps(intersect_types,ensure_ascii=False),'geometry':hole})
holes=gpd.GeoDataFrame(records,crs=32651)
# Filling would swallow another center: refuse instead of silently merging it.
for hole in holes.itertuples():
 for i in commercial.sindex.query(hole.geometry,predicate='intersects'):
  other=commercial.iloc[i]
  if int(other.id)!=hole.id and hole.geometry.intersection(other.geometry).area>.001:
   raise ValueError(f'商圈 {hole.id} 的内环包含商圈 {other.id}，需要人工处理')
commercial['geometry']=commercial.geometry.map(lambda x:Polygon(x.exterior))
commercial['area_km2']=commercial.area/1e6
commercial['hole_fill_m2']=commercial.id.map(holes.groupby('id').area_m2.sum()).fillna(0.)
fill=holes.geometry.union_all()
special['geometry']=special.geometry.map(lambda x: make_valid(x.difference(fill)))
special=special[~special.geometry.is_empty].copy();special['area_km2']=special.area/1e6
import pandas as pd
partition=gpd.GeoDataFrame(pd.concat([commercial,special],ignore_index=True),crs=32651)
area_union=partition.geometry.union_all(); commercial_union=commercial.geometry.union_all()
metrics={'trade_area_count':len(commercial),'hole_count_filled':len(holes),'affected_trade_areas':int((commercial.hole_fill_m2>0).sum()),'area_filled_km2':float(holes.area_m2.sum()/1e6),'hole_area_by_previous_zone_km2':{},'multipolygon_count':int((commercial.geom_type!='Polygon').sum()),'remaining_interior_rings':int(sum(len(g.interiors) for g in commercial.geometry)),'invalid_commercial_count':int((~commercial.is_valid).sum()),'commercial_overlap_m2':float(commercial.area.sum()-commercial_union.area),'partition_invalid_count':int((~partition.is_valid).sum()),'partition_overlap_m2':float(partition.area.sum()-area_union.area),'coverage_gap_m2':float(boundary.difference(area_union).area),'commercial_area_km2':float(commercial.area.sum()/1e6)}
for item in holes.contents:
 for item2 in json.loads(item):metrics['hole_area_by_previous_zone_km2'][item2['zone_type']]=metrics['hole_area_by_previous_zone_km2'].get(item2['zone_type'],0)+item2['area_m2']/1e6
assert metrics['multipolygon_count']==metrics['remaining_interior_rings']==metrics['invalid_commercial_count']==metrics['partition_invalid_count']==0
assert abs(metrics['commercial_overlap_m2'])<.1 and abs(metrics['partition_overlap_m2'])<.1 and metrics['coverage_gap_m2']<.1
# Coordinate reprojection can collapse tiny slivers into lines; retain only polygonal parts.
from shapely.ops import unary_union
def polygon_only(geom):
 geom=make_valid(geom)
 def parts(g):
  if g.geom_type=='Polygon': return [g]
  if hasattr(g,'geoms'): return [p for part in g.geoms for p in parts(part)]
  return []
 return unary_union(parts(geom))
commercial_wgs=commercial.to_crs(4326);partition_wgs=partition.to_crs(4326)
partition_wgs['geometry']=partition_wgs.geometry.map(polygon_only)
assert commercial_wgs.is_valid.all() and partition_wgs.is_valid.all()
commercial_wgs.to_file(out/'trade_areas_no_holes.geojson',driver='GeoJSON');partition_wgs.to_file(out/'city_partition_no_holes.geojson',driver='GeoJSON');holes.to_crs(4326).to_file(out/'filled_holes.geojson',driver='GeoJSON');(out/'metrics.json').write_text(json.dumps(metrics,indent=2,ensure_ascii=False));print(json.dumps(metrics,ensure_ascii=False,indent=2))
