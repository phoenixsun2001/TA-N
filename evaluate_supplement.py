from pathlib import Path
import json
import numpy as np,geopandas as g,shapely
from shapely.strtree import STRtree
out=Path('data/user_shanghai');city=g.read_file('data/raw/shanghai_boundary.geojson').to_crs(32651).geometry.union_all();shapely.prepare(city)
classes=['motorway','trunk','primary','secondary','tertiary','residential','unclassified','service','living_street']
def load(p):
 d=g.read_file(p).to_crs(32651);d=d[d.fclass.isin(classes)&(d.bridge!='T')&(d.tunnel!='T')&(d.layer.fillna(0)==0)];return d[shapely.intersects(city,d.geometry.array)].reset_index(drop=True)
new=load(out/'national_roads.gpkg');old=load('data/raw/osm/gis_osm_roads_free_1.shp');tree=STRtree(old.geometry.to_numpy());geom=new.geometry.to_numpy();mid=shapely.line_interpolate_point(geom,.5,normalized=True);_,d=tree.query_nearest(mid,return_distance=True,all_matches=False);candidate=np.flatnonzero((d>10)&(shapely.length(geom)>=10));print('midpoint candidates',len(candidate),flush=True)
selected=[];details=[]
for k in candidate:
 line=geom[k];n=max(3,int(np.ceil(line.length/20)));pts=shapely.line_interpolate_point(line,(np.arange(n)+.5)*line.length/n);pairs=tree.query(pts,predicate='dwithin',distance=5);near=len(np.unique(pairs[0]))/n
 if near<=.2:
  selected.append(int(k));details.append((float(d[k]),float(near)))
res=new.iloc[selected].copy();res['mid_dist']=np.array(details)[:,0] if details else [];res['near_share']=np.array(details)[:,1] if details else [];res['source']='national_2025_candidate';res.to_file(out/'national_roads_candidate_additions.gpkg',driver='GPKG');m={'source_national_ground_roads':len(new),'existing_ground_roads':len(old),'midpoint_candidates_over_10m':len(candidate),'selected_additions':len(res),'selected_length_km':float(res.length.sum()/1000),'selection':'midpoint >10 m, length >=10m, <=20% sampled points within 5m of existing road; not independently validated'};(out/'supplement_metrics.json').write_text(json.dumps(m,ensure_ascii=False,indent=2));print(json.dumps(m,indent=2),flush=True)
# Sample long hydrology lines only inside the Shanghai boundary / small buffer.
ref=g.read_file('data/raw/osm/gis_osm_waterways_free_1.shp').to_crs(32651);pol=g.read_file('data/raw/osm/gis_osm_water_a_free_1.shp').to_crs(32651);tree2=STRtree(np.concatenate([ref.geometry.to_numpy(),pol.geometry.to_numpy()]));buf=city.buffer(2000)
for fp in sorted(out.glob('我国*河流.shp.gpkg')):
 v=g.read_file(fp).to_crs(32651);seg=shapely.intersection(v.geometry.array,buf);length=shapely.length(seg);idx=np.flatnonzero(length>50);samples=[]
 for i in idx:
  n=max(1,int(np.ceil(length[i]/250)));samples.extend(shapely.line_interpolate_point(seg[i],(np.arange(n)+.5)*length[i]/n))
 if not samples:continue
 _,dist=tree2.query_nearest(np.array(samples,dtype=object),return_distance=True,all_matches=False);print(fp.name,'in-range samples',len(samples),'median',round(float(np.median(dist)),1),'p90',round(float(np.quantile(dist,.9)),1),flush=True)
