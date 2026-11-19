from pathlib import Path
import json,time,warnings
import geopandas as gpd,numpy as np,shapely
from shapely.strtree import STRtree
from convert_gcj02 import convert
out=Path('outputs/rebuild_v6');out.mkdir(parents=True,exist_ok=True)
base=None
city=gpd.read_file('data/raw/shanghai_boundary.geojson').to_crs(32651).geometry.union_all()
base=gpd.GeoDataFrame({'unit_seed':[0]},geometry=[city],crs=32651)
shapely.prepare(city)
# Estimated half-widths in metres; centerlines do not encode actual right-of-way.
half={'motorway':15,'motorway_link':8,'trunk':12,'trunk_link':7,'primary':9,'primary_link':6,'secondary':7,'secondary_link':5,'tertiary':6,'tertiary_link':5,'residential':4,'unclassified':4,'living_street':3,'service':3,'track':2,'track_grade1':2,'track_grade2':2,'track_grade3':2,'track_grade5':2}
roads=gpd.read_file('data/raw/osm/gis_osm_roads_free_1.shp').to_crs(32651)
extra=gpd.read_file('data/user_shanghai/national_roads_candidate_filtered.gpkg').to_crs(32651)
import pandas as pd
roads=gpd.GeoDataFrame(pd.concat([roads,extra],ignore_index=True),crs=32651)
selected=['motorway','motorway_link','trunk','trunk_link','primary','primary_link','secondary','secondary_link','tertiary','tertiary_link']
roads=roads[roads.fclass.isin(selected)&(roads.tunnel!='T')&(roads.layer.fillna(0)>=0)].copy()
roads=roads[shapely.intersects(city,roads.geometry.array)].copy(); roads['half_width_m']=roads.fclass.map(half)
roadbuf=shapely.buffer(roads.geometry.to_numpy(),roads.half_width_m.to_numpy(),quad_segs=2)
rail=gpd.read_file('data/raw/osm/gis_osm_railways_free_1.shp').to_crs(32651)
rail=rail[(rail.fclass.isin(['rail','light_rail','tram','monorail','narrow_gauge','funicular'])|((rail.fclass=='subway')&((rail.bridge=='T')|(rail.layer.fillna(0)>0))))&(rail.tunnel!='T')&(rail.layer.fillna(0)>=0)].copy()
rail=rail[shapely.intersects(city,rail.geometry.array)].copy();railbuf=shapely.buffer(rail.geometry.to_numpy(),5,quad_segs=2)
water=gpd.read_file('data/raw/osm/gis_osm_water_a_free_1.shp').to_crs(32651)
water=water[water.fclass.isin(['water','riverbank','reservoir','dock'])&shapely.intersects(city,water.geometry.array)]
watergeom=shapely.make_valid(water.geometry.to_numpy())
trees=[STRtree(roadbuf),STRtree(railbuf),STRtree(watergeom)]; arrays=[roadbuf,railbuf,watergeom]
print('inputs',len(base),len(roads),len(rail),len(water),flush=True)
rows=[];before=0.;after=0.;excluded=0.;fully_removed=0;start=time.time()
for idx,row in base.iterrows():
 geom=row.geometry;before+=geom.area
 near=[a[t.query(geom,predicate='intersects')] for a,t in zip(arrays,trees)]
 mask=shapely.union_all(np.concatenate(near)) if any(len(a) for a in near) else None
 cut=geom.difference(mask) if mask is not None else geom
 if not cut.is_valid:cut=shapely.make_valid(cut)
 stack=[cut];parts=[]
 while stack:
  z=stack.pop()
  if z.geom_type=='Polygon' and not z.is_empty and z.area>0:parts.append(z)
  elif hasattr(z,'geoms'):stack.extend(z.geoms)
 if not parts:fully_removed+=1
 for k,p in enumerate(parts):
  attrs=row.drop(labels='geometry').to_dict();attrs['parent_seed']=int(attrs.pop('unit_seed'));attrs['component']=k+1;attrs['area_km2']=p.area/1e6;attrs['geometry']=p;rows.append(attrs);after+=p.area
 if (idx+1)%1000==0:print('processed',idx+1,'components',len(rows),'seconds',round(time.time()-start),flush=True)
res=gpd.GeoDataFrame(rows,crs=32651);res.insert(0,'block_id',np.arange(1,len(res)+1,dtype='int32'))
res.to_file(out/'land_blocks_projected.gpkg',driver='GPKG')
print('saving GCJ',len(res),flush=True)
gcj=res.to_crs(4326);gcj.geometry=shapely.transform(gcj.geometry.array,convert,include_z=False);gcj.crs=None
# Nonlinear conversion can cause tiny self-intersections. Repair and split the
# resulting polygonal parts; retain block_id as the source projected ID.
bad=~gcj.is_valid
gcj.loc[bad,'geometry']=shapely.make_valid(gcj.loc[bad,'geometry'].array)
gcj=gcj.explode(index_parts=False,ignore_index=True)
gcj=gcj[gcj.geom_type=='Polygon'].copy()
gcj['gcj_part']=gcj.groupby('block_id').cumcount()+1
gcj.insert(0,'feature_id',np.arange(1,len(gcj)+1,dtype='int32'))
gcj=gcj.rename(columns={'area_km2':'source_area_km2'})
with warnings.catch_warnings():
 warnings.filterwarnings('ignore',message="'crs' was not provided.*");gcj.to_file(out/'land_blocks_gcj02.geojson',driver='GeoJSON')
# Inspect overlap and interior residual inside masks; excludes boundaries by area.
invalid=int((~res.is_valid).sum());holes=int(sum(len(z.interiors)>0 for z in res.geometry));multis=int((res.geom_type!='Polygon').sum())
metrics={'input_units':len(base),'output_components':len(res),'source_roads':len(roads),'source_rail_lines':len(rail),'source_water_polygons':len(water),'estimated_area_before_km2':before/1e6,'area_after_km2':after/1e6,'excluded_area_km2':(before-after)/1e6,'fully_removed_input_units':fully_removed,'invalid':invalid,'non_polygon_components':multis,'components_with_holes':holes,'road_half_width_m':{k:half[k] for k in selected},'rail_half_width_m':5,'note':'Cuts motorway through tertiary and their links, including elevated roads; excludes tunnels and negative layers. Lower road classes are not removed.'}
(out/'metrics.json').write_text(json.dumps(metrics,ensure_ascii=False,indent=2));print(json.dumps(metrics,indent=2,ensure_ascii=False),flush=True)
