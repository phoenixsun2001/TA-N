from pathlib import Path
import json
import numpy as np,geopandas as gpd,pandas as pd,shapely
from shapely.geometry import Polygon
from shapely.ops import unary_union
from shapely.strtree import STRtree
out=Path('outputs/block_fit');b=gpd.read_file(out/'city_block_assignment.gpkg').to_crs(32651);anchors=gpd.read_file(out/'city_seed_anchors.gpkg').to_crs(32651);city=gpd.read_file('data/raw/shanghai_boundary.geojson').to_crs(32651).geometry.union_all();assigned=b[b.seed_id>0];parts=[];holes=0
for seed,group in assigned.groupby('seed_id'):
 original=unary_union(group.geometry);p=list(original.geoms) if original.geom_type=='MultiPolygon' else [original]
 holes+=sum(len(g.interiors) for g in p if g.geom_type=='Polygon')
 filled=unary_union([Polygon(g.exterior) for g in p if g.geom_type=='Polygon'])
 for k,poly in enumerate(filled.geoms if filled.geom_type=='MultiPolygon' else [filled],1):
  parts.append({'area_id':f'{seed}-{k:02d}','seed_id':int(seed),'geometry':poly})
trade=gpd.GeoDataFrame(parts,crs=32651);tree=STRtree(trade.geometry.to_numpy());pairs=tree.query(trade.geometry.to_numpy(),predicate='intersects').T
conflicts=[(int(i),int(j)) for i,j in pairs if i<j and trade.geometry.iloc[i].intersection(trade.geometry.iloc[j]).area>.01]
if conflicts:raise RuntimeError(f'填充内环后商圈互相重叠，不能自动输出: {conflicts[:5]}')
trade['area_km2']=trade.area/1e6;trade_union=trade.geometry.union_all()
# Keep water overlay and unassigned complementary so the full input city remains covered.
water=b[b.is_water].geometry.union_all().intersection(city).difference(trade_union);unassigned=city.difference(unary_union([water,trade_union]));special=gpd.GeoDataFrame([{'area_id':'water','seed_id':0,'zone_type':'water','geometry':water},{'area_id':'unassigned','seed_id':0,'zone_type':'unassigned','geometry':unassigned}],crs=32651);trade['zone_type']='commercial_candidate';part=gpd.GeoDataFrame(pd.concat([trade,special],ignore_index=True),crs=32651)
def polygon_only(geom):
 def parts(g):
  if g.geom_type=='Polygon':return [g]
  if hasattr(g,'geoms'):return [p for x in g.geoms for p in parts(x)]
  return []
 return unary_union(parts(shapely.make_valid(geom)))
wgs=part.to_crs(4326);wgs.geometry=wgs.geometry.map(polygon_only);tw=trade.to_crs(4326)
assert tw.is_valid.all() and (tw.geom_type=='Polygon').all() and sum(len(g.interiors) for g in tw.geometry)==0 and wgs.is_valid.all()
tw.to_file(out/'city_street_trade_areas.geojson',driver='GeoJSON');wgs.to_file(out/'city_street_partition.geojson',driver='GeoJSON');anchors.to_crs(4326).to_file(out/'city_seed_anchors.geojson',driver='GeoJSON')
# Candidate source roads; estimate actual distance of sampled output boundary to road or other mapped barrier.
roads=gpd.read_file('data/raw/osm/gis_osm_roads_free_1.shp').to_crs(32651);cut=roads[roads.fclass.isin(['motorway','trunk','primary','secondary','tertiary','residential','unclassified','service','living_street'])&(roads.bridge!='T')&(roads.tunnel!='T')&(roads.layer.fillna(0)==0)];rtree=STRtree(cut.geometry.to_numpy());major_tree=STRtree(cut[cut.fclass.isin(['motorway','trunk','primary','secondary'])].geometry.to_numpy())
near=near_major=tot=0.
for geom in trade.geometry:
 line=geom.boundary;length=line.length;n=max(1,int(np.ceil(length/25)));pts=shapely.line_interpolate_point(line,(np.arange(n)+.5)*length/n);hits=rtree.query(pts,predicate='dwithin',distance=3);near+=len(np.unique(hits[0]))*length/n;majorhits=major_tree.query(pts,predicate='dwithin',distance=3);near_major+=len(np.unique(majorhits[0]))*length/n;tot+=length
check=gpd.read_file(out/'city_street_partition.geojson').to_crs(32651);u=check.geometry.union_all();cc=check[check.seed_id>0];missing=[]
for r in anchors[anchors.matched].itertuples():
 p=cc[cc.seed_id==r.seed_id].geometry.union_all()
 if p.distance(r.geometry)>0.05:missing.append(int(r.seed_id))
metrics=json.loads((out/'city_graph_metrics.json').read_text());metrics.update({'components':len(trade),'original_interior_rings_filled':holes,'commercial_area_filled_km2':float(trade.area.sum()/1e6),'road_boundary_within_3m_pct':near/tot*100,'major_road_boundary_within_3m_pct':near_major/tot*100,'invalid_partition':int((~check.is_valid).sum()),'commercial_holes':int(sum(len(g.interiors) for g in cc.geometry)),'commercial_nonpolygons':int((cc.geom_type!='Polygon').sum()),'commercial_overlap_m2':float(cc.area.sum()-cc.geometry.union_all().area),'partition_overlap_m2':float(check.area.sum()-u.area),'partition_gap_m2':float(city.difference(u).area),'anchors_outside_parent_geometry':missing})
assert metrics['invalid_partition']==metrics['commercial_holes']==metrics['commercial_nonpolygons']==0 and abs(metrics['commercial_overlap_m2'])<.1 and abs(metrics['partition_overlap_m2'])<.1 and metrics['partition_gap_m2']<.1 and len(missing)==0
(out/'city_final_metrics.json').write_text(json.dumps(metrics,ensure_ascii=False,indent=2));print(json.dumps(metrics,ensure_ascii=False,indent=2),flush=True)
