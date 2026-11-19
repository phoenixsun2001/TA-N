import json
from pathlib import Path
import numpy as np,geopandas as gpd,shapely
from shapely.strtree import STRtree
p=Path('outputs/barrier_cut_v3');d=gpd.read_file(p/'land_blocks_projected.gpkg');n=len(d);areas=d.area.to_numpy();print('read',n,flush=True)
road=gpd.read_file('data/raw/osm/gis_osm_roads_free_1.shp').to_crs(32651);extra=gpd.read_file('data/user_shanghai/national_roads_candidate_filtered.gpkg').to_crs(32651)
import pandas as pd
road=gpd.GeoDataFrame(pd.concat([road,extra],ignore_index=True),crs=32651)
half={'motorway':15,'motorway_link':8,'trunk':12,'trunk_link':7,'primary':9,'primary_link':6,'secondary':7,'secondary_link':5,'tertiary':6,'tertiary_link':5,'residential':4,'unclassified':4,'living_street':3,'service':3,'track':2,'track_grade1':2,'track_grade2':2,'track_grade3':2,'track_grade5':2}
road=road[road.fclass.isin(half)&(road.tunnel!='T')&(road.layer.fillna(0)>=0)];rg=shapely.buffer(road.geometry.to_numpy(),road.fclass.map(half).to_numpy(),quad_segs=2)
rail=gpd.read_file('data/raw/osm/gis_osm_railways_free_1.shp').to_crs(32651);rail=rail[(rail.fclass.isin(['rail','light_rail','tram','monorail','narrow_gauge','funicular'])|((rail.fclass=='subway')&((rail.bridge=='T')|(rail.layer.fillna(0)>0))))&(rail.tunnel!='T')&(rail.layer.fillna(0)>=0)];rr=shapely.buffer(rail.geometry.to_numpy(),5,quad_segs=2)
water=gpd.read_file('data/raw/osm/gis_osm_water_a_free_1.shp').to_crs(32651);water=water[water.fclass.isin(['water','riverbank','reservoir','dock'])];ww=shapely.make_valid(water.geometry.array)
trees=[STRtree(a) for a in [rg,rr,ww]];arrays=[rg,rr,ww];rng=np.random.default_rng(12345);sample=rng.choice(n,min(n,1000),replace=False);residual=[0.,0.,0.]
for i in sample:
 z=d.geometry.iloc[i]
 for j,(tree,arr) in enumerate(zip(trees,arrays)):
  a=arr[tree.query(z,predicate='intersects')]
  if len(a):residual[j]+=z.intersection(shapely.union_all(a)).area
print('sample test done',flush=True)
gcj=gpd.read_file(p/'land_blocks_gcj02.geojson');print('gcj read',flush=True)
m={'components':n,'small_under_100_m2':int((areas<100).sum()),'small_under_1000_m2':int((areas<1000).sum()),'area_under_1000_m2_km2':float(areas[areas<1000].sum()/1e6),'zone_types':d.zone_type.value_counts().to_dict(),'unresolved_component_area_km2':float(d.area[d.zone_type=='road_network_incomplete'].sum()/1e6),'sample_count':len(sample),'sample_residual_overlap_m2':dict(zip(['road','rail','water'],residual)),'gcj_features':len(gcj),'gcj_invalid':int((~gcj.is_valid).sum()),'gcj_holes':int(sum(z.geom_type=='Polygon' and len(z.interiors)>0 for z in gcj.geometry))}
(p/'verification.json').write_text(json.dumps(m,indent=2));print(json.dumps(m,indent=2),flush=True)
