from pathlib import Path
import json,heapq,time
import numpy as np,geopandas as gpd,shapely,rasterio
from rasterio.features import rasterize
from shapely.geometry import Point
from shapely.strtree import STRtree
from scipy.ndimage import median_filter
out=Path('outputs/block_fit')
def read(name):return gpd.read_file(f'data/raw/osm/gis_osm_{name}_free_1.shp').to_crs(32651)
blocks=gpd.read_file(out/'city_blocks.gpkg').to_crs(32651);geom=blocks.geometry.to_numpy();n=len(blocks);land=~blocks.is_water.to_numpy();big=blocks.area.to_numpy()>2e6
roads=read('roads');rails=read('railways');traffic=read('traffic');centers=gpd.read_file('outputs/centers.geojson').to_crs(32651)
major=roads[roads.fclass.isin(['motorway','trunk','primary','secondary'])&(roads.bridge!='T')&(roads.tunnel!='T')&(roads.layer.fillna(0)==0)];rail=rails[rails.fclass.isin(['rail','light_rail','tram'])&(rails.bridge!='T')&(rails.tunnel!='T')&(rails.layer.fillna(0)==0)];walk=roads[roads.fclass.isin(['footway','path','pedestrian','steps'])];bridge=walk[(walk.bridge=='T')|(walk.tunnel=='T')];ped=traffic[traffic.fclass=='pedestrian_crossing'];level=traffic[traffic.fclass=='railway_crossing']
majortree=STRtree(major.geometry.to_numpy());railtree=STRtree(rail.geometry.to_numpy());wtree=STRtree(walk.geometry.to_numpy());ptree=STRtree(ped.geometry.to_numpy());ltree=STRtree(level.geometry.to_numpy());btree=STRtree(geom)
print('blocks',n,'land',land.sum(),'large',big.sum(),flush=True)
# Attribute available satellite context on the 250 m feature grid, without treating resampled VIIRS as independent detail.
with rasterio.open('outputs/city_imagery/city_features.tif') as ds:
 cell=rasterize(((g,i+1) for i,g in enumerate(geom)),out_shape=(ds.height,ds.width),transform=ds.transform,fill=0,dtype='int32')-1
 features={}
 for name in ['viirs_index','ndvi','ndbi','s1_vv_db']:
  values=ds.read(ds.descriptions.index(name)+1);ok=(cell>=0)&(values!=-9999)&np.isfinite(values)
  count=np.bincount(cell[ok],minlength=n);sums=np.bincount(cell[ok],weights=values[ok],minlength=n);features[name]=np.divide(sums,count,out=np.full(n,np.nan),where=count>0)
blocks['ntl']=features['viirs_index'];blocks['ndvi']=features['ndvi'];blocks['ndbi']=features['ndbi'];blocks['s1_vv']=features['s1_vv_db']
ntl=features['viirs_index'];ndvi=features['ndvi'];ndbi=features['ndbi'];vv=features['s1_vv_db'];factor=np.where(np.isfinite(ntl)&(ntl<3),1.15,np.where(np.isfinite(ntl)&(ntl>21),.90,1.0))*np.where(np.isfinite(ndvi)&(ndvi>.35),1.8,1.0)*np.where(np.isfinite(vv)&(vv>-8)&np.isfinite(ndbi)&(ndbi>0),.90,1.0)
# Shared-edge street-block graph; crossing indicators are proxies, not a verified walk network.
reps=np.array([g.representative_point() for g in geom]);pairs=btree.query(geom,predicate='intersects').T;adj=[[] for _ in range(n)];records=[];counts={'ordinary':0,'major_with_crossing':0,'rail_with_crossing':0,'blocked_major':0,'blocked_rail':0,'blocked_large':0,'bridge_link':0}
def linked(i,j,cost,kind,geometry,source=''):
 adj[i].append((j,cost));adj[j].append((i,cost));records.append({'block_a':i,'block_b':j,'kind':kind,'cost_m_proxy':cost,'source_osm_id':source,'geometry':geometry})
def covered_fraction(shared,tree,geo):
 ix=tree.query(shared,predicate='dwithin',distance=3)
 if len(ix)==0:return 0.0
 near=shapely.union_all(geo[ix]);return shared.intersection(near.buffer(3)).length/shared.length
mg=major.geometry.to_numpy();rg=rail.geometry.to_numpy();wg=walk.geometry.to_numpy()
for i,j in pairs:
 i=int(i);j=int(j)
 if i>=j or not land[i] or not land[j]:continue
 if big[i] or big[j]:counts['blocked_large']+=1;continue
 shared=geom[i].boundary.intersection(geom[j].boundary)
 if shared.length<1:continue
 onrail=covered_fraction(shared,railtree,rg)>=.3;onmajor=covered_fraction(shared,majortree,mg)>=.3
 level_cross=len(ltree.query(shared,predicate='dwithin',distance=20))>0 if onrail else False
 ped_cross=len(ptree.query(shared,predicate='dwithin',distance=20))>0 if onmajor else False
 walking_cross=any(wg[k].crosses(shared) for k in wtree.query(shared,predicate='intersects')) if onrail or onmajor else False
 if onrail and not(level_cross or walking_cross):counts['blocked_rail']+=1;continue
 if onmajor and not(ped_cross or walking_cross):counts['blocked_major']+=1;continue
 kind='rail_with_crossing' if onrail else 'major_with_crossing' if onmajor else 'ordinary';counts[kind]+=1
 linked(i,j,float(reps[i].distance(reps[j]))+(100 if onrail or onmajor else 30),kind,shared)
print('planar',counts,flush=True)
water_tree=STRtree(geom[~land]);watergeom=geom[~land]
for way in bridge.itertuples():
 if way.geometry.geom_type!='LineString' or way.geometry.length<10:continue
 # Keep only explicitly tagged foot bridges/tunnels that actually intersect a mapped barrier.
 barrier=bool(len(majortree.query(way.geometry,predicate='intersects')) or len(railtree.query(way.geometry,predicate='intersects')) or len(water_tree.query(way.geometry,predicate='intersects')))
 if not barrier:continue
 sides=[];pts=[Point(way.geometry.coords[0]),Point(way.geometry.coords[-1])]
 for p in pts:
  ix=[int(i) for i in btree.query(p,predicate='dwithin',distance=35) if land[int(i)] and not big[int(i)]]
  if not ix:sides=[];break
  sides.append(min(ix,key=lambda i:geom[i].distance(p)))
 if len(sides)!=2 or sides[0]==sides[1]:continue
 i,j=sides
 if geom[i].boundary.intersection(geom[j].boundary).length>=1:continue
 linked(i,j,float(reps[i].distance(pts[0])+way.geometry.length+reps[j].distance(pts[1])+100),'bridge_link',way.geometry,str(way.osm_id));counts['bridge_link']+=1
print('links',counts,flush=True)
seednodes=[];collisions=[];miss=[];used=set();anchor_records=[]
poi=gpd.read_file('outputs/commercial_poi.geojson').to_crs(32651);poigeom=poi.geometry.to_numpy();poitree=STRtree(poigeom)
for r in centers.itertuples():
 anchor=r.geometry;source='grid_center';displacement=0.0;source_poi=''
 candidates=[int(i) for i in btree.query(anchor,predicate='intersects') if land[int(i)] and not big[int(i)] and int(i) not in used]
 if not candidates:
  nearby=poitree.query(r.geometry,predicate='dwithin',distance=150);nearby=sorted(nearby,key=lambda i:poigeom[i].distance(r.geometry))
  for k in nearby:
   pt=poigeom[k];found=[int(i) for i in btree.query(pt,predicate='intersects') if land[int(i)] and not big[int(i)] and int(i) not in used]
   if found:
    candidates=found;anchor=pt;source='nearby_commercial_poi';displacement=float(pt.distance(r.geometry));source_poi=str(poi.osm_id.iloc[k]);break
 if not candidates:
  miss.append(int(r.id));anchor_records.append({'seed_id':int(r.id),'matched':False,'source':'unmatched','distance_m':-1.0,'source_poi':'','geometry':r.geometry});continue
 node=min(candidates,key=lambda i:geom[i].distance(anchor));used.add(node);seednodes.append((int(r.id),node))
 anchor_records.append({'seed_id':int(r.id),'matched':True,'source':source,'distance_m':displacement,'source_poi':source_poi,'geometry':anchor})
gpd.GeoDataFrame(anchor_records,crs=32651).to_file(out/'city_seed_anchors.gpkg',driver='GPKG')
def solve(costfactor):
 dist=np.full(n,np.inf);labels=np.zeros(n,dtype='int32');heap=[]
 for id,node in seednodes:dist[node]=0;labels[node]=id;heapq.heappush(heap,(0.,id,node))
 while heap:
  d,id,i=heapq.heappop(heap)
  if d>dist[i]+1e-7 or labels[i]!=id:continue
  for j,c in adj[i]:
   nd=d+c*(costfactor[i]+costfactor[j])/2
   if nd<=2000 and nd<dist[j]:dist[j]=nd;labels[j]=id;heapq.heappush(heap,(nd,id,j))
 return labels,dist
base,_=solve(np.ones(n));labels,dist=solve(factor);blocks['seed_id']=labels;blocks['baseline_id']=base;blocks['cost_proxy']=np.where(np.isfinite(dist),dist,-1);blocks['large_flag']=big;blocks.to_file(out/'city_block_assignment.gpkg',driver='GPKG')
if records:gpd.GeoDataFrame(records,crs=32651).to_file(out/'city_passage_links.gpkg',driver='GPKG')
result={'blocks':n,'land_blocks':int(land.sum()),'large_blocks':int(big.sum()),'seed_total':len(centers),'seeds_matched':len(seednodes),'seeds_missing':miss,'seeds_colliding':collisions,'seeds_relocated_to_poi':sum(r['source']=='nearby_commercial_poi' for r in anchor_records),'assigned_blocks':int((labels>0).sum()),'unassigned_land_blocks':int(((labels==0)&land).sum()),'assigned_area_km2':float(blocks.area[labels>0].sum()/1e6),'changed_blocks':int((base!=labels).sum()),'graph_edges':counts,'blocks_with_s2_context':int(np.isfinite(ndvi).sum()),'blocks_with_viirs_context':int(np.isfinite(ntl).sum())}
(out/'city_graph_metrics.json').write_text(json.dumps(result,ensure_ascii=False,indent=2));print(json.dumps(result,ensure_ascii=False,indent=2),flush=True)
