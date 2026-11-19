from pathlib import Path
import json,heapq,time
import pandas as pd
import numpy as np,geopandas as gpd,shapely,rasterio
from rasterio.features import rasterize
from shapely.geometry import Point
from shapely.strtree import STRtree
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
out=Path('outputs/user_data_test')
def read(name):return gpd.read_file(f'data/raw/osm/gis_osm_{name}_free_1.shp').to_crs(32651)
blocks=gpd.read_file(out/'candidate_blocks.gpkg').to_crs(32651);geom=blocks.geometry.to_numpy();n=len(blocks);land=~blocks.is_water.to_numpy();big=blocks.area.to_numpy()>2e6
roads=read('roads');rails=read('railways');traffic=read('traffic');centers=gpd.read_file('outputs/centers.geojson').to_crs(32651)
supplement=gpd.read_file('data/user_shanghai/national_roads_candidate_filtered.gpkg').to_crs(32651);major=gpd.GeoDataFrame(pd.concat([roads[roads.fclass.isin(['motorway','trunk','primary','secondary'])&(roads.bridge!='T')&(roads.tunnel!='T')&(roads.layer.fillna(0)==0)],supplement[supplement.fclass.isin(['motorway','trunk','primary','secondary'])]],ignore_index=True),crs=32651);rail=rails[rails.fclass.isin(['rail','light_rail','tram'])&(rails.bridge!='T')&(rails.tunnel!='T')&(rails.layer.fillna(0)==0)];walk=roads[roads.fclass.isin(['footway','path','pedestrian','steps'])];bridge=walk[(walk.bridge=='T')|(walk.tunnel=='T')];ped=traffic[traffic.fclass=='pedestrian_crossing'];level=traffic[traffic.fclass=='railway_crossing']
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
reps=np.array([g.representative_point() for g in geom]);eligible_ids=np.flatnonzero(land&~big);eligible_tree=STRtree(geom[eligible_ids]);raw_pairs=eligible_tree.query(geom[eligible_ids],predicate='intersects').T;pairs=eligible_ids[raw_pairs];adj=[[] for _ in range(n)];records=[];counts={'ordinary':0,'major_with_crossing':0,'rail_with_crossing':0,'blocked_major':0,'blocked_rail':0,'blocked_large':0,'bridge_link':0,'inferred_unverified':0}
def linked(i,j,cost,kind,geometry,source=''):
 adj[i].append((j,cost));adj[j].append((i,cost));records.append({'block_a':i,'block_b':j,'kind':kind,'cost_m_proxy':cost,'source_osm_id':source,'geometry':geometry})
def covered_fraction(shared,tree,geo):
 ix=tree.query(shared,predicate='dwithin',distance=3)
 if len(ix)==0:return 0.0
 near=shapely.union_all(geo[ix]);return shared.intersection(near.buffer(3)).length/shared.length
mg=major.geometry.to_numpy();rg=rail.geometry.to_numpy();wg=walk.geometry.to_numpy();free=major[major.fclass.isin(['motorway','trunk'])].geometry.to_numpy();ftree=STRtree(free)
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
 inferred=onmajor and not(ped_cross or walking_cross)
 if inferred and covered_fraction(shared,ftree,free)>=.3:counts['blocked_major']+=1;continue
 kind='inferred_unverified' if inferred else 'rail_with_crossing' if onrail else 'major_with_crossing' if onmajor else 'ordinary';counts[kind]+=1
 linked(i,j,float(reps[i].distance(reps[j]))+(600 if inferred else 100 if onrail or onmajor else 30),kind,shared)
counts['blocked_large']=None;print('planar',counts,flush=True)
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
# Commercial AOI/POI anchors are independent of full-city activity anchors.
seed_by_node={};seed_info=[];mall=gpd.read_file('data/raw/osm/gis_osm_pois_a_free_1.shp').to_crs(32651);city=gpd.read_file('data/raw/shanghai_boundary.geojson').to_crs(32651).geometry.union_all();mall=mall[(mall.fclass=='mall')&mall.geometry.representative_point().within(city)]
def locate(p):
 candidates=[int(i) for i in btree.query(p,predicate='intersects') if land[int(i)] and not big[int(i)] and geom[int(i)].covers(p)]
 return min(candidates,key=lambda i:geom[i].area) if candidates else None
def put(node,kind,label):
 if node not in seed_by_node:
  sid=len(seed_info)+1;seed_by_node[node]=sid;seed_info.append({'seed_id':sid,'seed_kind':kind,'anchor_labels':[],'seed_block':node})
 sid=seed_by_node[node];seed_info[sid-1]['anchor_labels'].append(str(label));return sid
mall_status=[]
for row in mall.itertuples():
 point=row.geometry.representative_point();node=locate(point)
 if node is not None:put(node,'commercial_aoi',row.name or row.osm_id)
 mall_status.append({'osm_id':str(row.osm_id),'name':row.name,'matched':node is not None,'geometry':point})
for row in centers.itertuples():
 node=locate(row.geometry)
 if node is not None:put(node,'poi_density',row.id)
supp=gpd.read_file('outputs/revised_poc/additional_poi_clusters.geojson').to_crs(32651)
for row in supp.itertuples():
 for part in shapely.get_parts(row.geometry):
  if part.geom_type!='Polygon':continue
  node=locate(part.representative_point())
  if node is not None:put(node,'poi_cluster',row.source_id)
commercial_count=len(seed_info);print('commercial anchors',commercial_count,'mall matched',sum(r['matched'] for r in mall_status),'/',len(mall_status),flush=True)
active=land&~big;seed_bins={(int(np.floor(reps[node].x/1500)),int(np.floor(reps[node].y/1500))) for node in seed_by_node};grid={}
for node in np.flatnonzero(active):
 p=reps[node];key=(int(np.floor(p.x/1500)),int(np.floor(p.y/1500)))
 if key in seed_bins:continue
 d=(p.x-(key[0]+.5)*1500)**2+(p.y-(key[1]+.5)*1500)**2
 if key not in grid or d<grid[key][0]:grid[key]=(d,int(node))
for _,node in grid.values():put(node,'local','')
rr=[];cc=[]
for i,edges in enumerate(adj):
 for j,_ in edges:rr.append(i);cc.append(j)
mat=coo_matrix((np.ones(len(rr),dtype='int8'),(rr,cc)),shape=(n,n)).tocsr();_,comp=connected_components(mat,directed=False);active_nodes=np.flatnonzero(active);anchored=set(comp[list(seed_by_node)]);missing=set(comp[active_nodes])-anchored;order=np.argsort(comp[active_nodes],kind='stable');sorted_nodes=active_nodes[order];cs=comp[sorted_nodes];starts=np.r_[0,np.flatnonzero(np.diff(cs))+1];ends=np.r_[starts[1:],len(sorted_nodes)];areas=blocks.area.to_numpy()
for start,end in zip(starts,ends):
 if int(cs[start]) not in missing:continue
 nodes=sorted_nodes[start:end];put(int(nodes[np.argmax(areas[nodes])]),'local_component','')
print('seeds',len(seed_info),'extra disconnected components',len(missing),flush=True)
dist=np.full(n,np.inf);labels=np.zeros(n,dtype='int32');queue=[]
for node,sid in seed_by_node.items():dist[node]=0;labels[node]=sid;heapq.heappush(queue,(0.,sid,node))
while queue:
 d,sid,i=heapq.heappop(queue)
 if d>dist[i]+1e-7 or labels[i]!=sid:continue
 for j,cost in adj[i]:
  nd=d+cost*(factor[i]+factor[j])/2
  if nd<dist[j]:dist[j]=nd;labels[j]=sid;heapq.heappush(queue,(nd,sid,j))
assert np.all(labels[active]>0)
blocks['unit_seed']=labels;blocks['cost_proxy']=np.where(np.isfinite(dist),dist,-1);blocks['large_flag']=big;blocks.to_file(out/'repaired_block_assignment.gpkg',driver='GPKG')
if records:gpd.GeoDataFrame(records,crs=32651).to_file(out/'repaired_passage_links.gpkg',driver='GPKG')
gpd.GeoDataFrame(mall_status,crs=32651).to_crs(4326).to_file(out/'mall_anchor_audit.geojson',driver='GeoJSON');(out/'road_seed_info.json').write_text(json.dumps(seed_info,ensure_ascii=False))
metrics={'blocks':n,'eligible_land_faces':int(active.sum()),'unresolved_large_faces':int((land&big).sum()),'matched_mall_anchors':sum(r['matched'] for r in mall_status),'total_mall_anchors':len(mall_status),'seed_count':len(seed_info),'commercial_seeds':commercial_count,'local_seeds':len(seed_info)-commercial_count,'unassigned_eligible_faces':int((active&(labels==0)).sum()),'graph_edges':counts};(out/'road_graph_metrics.json').write_text(json.dumps(metrics,ensure_ascii=False,indent=2));print(json.dumps(metrics,ensure_ascii=False,indent=2),flush=True)
