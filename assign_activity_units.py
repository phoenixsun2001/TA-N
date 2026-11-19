from pathlib import Path
import json,heapq
import numpy as np,geopandas as g,shapely,rasterio
from shapely.geometry import Point
from shapely.strtree import STRtree
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
out=Path('outputs/revised_poc');b=g.read_file(out/'activity_blocks.gpkg').to_crs(32651);geo=b.geometry.to_numpy();n=len(b);active=b.active.to_numpy()&~b.is_water.to_numpy();child=b.synthetic_split.to_numpy()&active;reps=np.array([x.representative_point() for x in geo]);active_ids=np.flatnonzero(active);idx=STRtree(geo[active_ids]);print('active land',int(active.sum()),'fallback',int(child.sum()),flush=True)
# Reuse already verified small-face links; retired oversized faces are filtered out.
links=g.read_file('outputs/block_fit/city_passage_links.gpkg').to_crs(32651);adj=[[] for _ in range(n)];edges=[]
def add(i,j,cost,kind):
 if i==j or not active[i] or not active[j]:return
 cost=max(float(cost),1.);adj[i].append((j,cost));adj[j].append((i,cost));edges.append((i,j,kind))
for r in links.itertuples():add(int(r.block_a),int(r.block_b),r.cost_m_proxy,r.kind)
print('reused edges',len(edges),flush=True)
def read(nm):return g.read_file(f'data/raw/osm/gis_osm_{nm}_free_1.shp').to_crs(32651)
roads=read('roads');rails=read('railways');traffic=read('traffic')
major=roads[roads.fclass.isin(['motorway','trunk','primary','secondary'])&(roads.bridge!='T')&(roads.tunnel!='T')&(roads.layer.fillna(0)==0)];rail=rails[rails.fclass.isin(['rail','light_rail','tram'])&(rails.bridge!='T')&(rails.tunnel!='T')&(rails.layer.fillna(0)==0)];walk=roads[roads.fclass.isin(['footway','path','pedestrian','steps'])];ped=traffic[traffic.fclass=='pedestrian_crossing'];level=traffic[traffic.fclass=='railway_crossing']
mg=major.geometry.to_numpy();rg=rail.geometry.to_numpy();wg=walk.geometry.to_numpy();mt=STRtree(mg);rt=STRtree(rg);wt=STRtree(wg);pt=STRtree(ped.geometry.to_numpy());lt=STRtree(level.geometry.to_numpy());free=major[major.fclass.isin(['motorway','trunk'])].geometry.to_numpy();ft=STRtree(free);inferred=[]
def frac(shared,tree,geoms):
 ids=tree.query(shared,predicate='dwithin',distance=3)
 return 0. if not len(ids) else shared.intersection(shapely.union_all(geoms[ids]).buffer(3)).length/shared.length
children=np.where(child)[0];pairs=idx.query(geo[children],predicate='intersects').T;seen=set();blocked={'major':0,'rail':0};newcount=0
for ci,j in pairs:
 i=int(children[ci]);j=int(active_ids[j])
 if not active[j] or i==j:continue
 key=(min(i,j),max(i,j))
 if key in seen:continue
 seen.add(key);shared=geo[i].boundary.intersection(geo[j].boundary)
 if shared.length<1:continue
 # Deliberate fallback grid edges inside the same large face are traversable.
 if child[i] and child[j] and b.old_id.iloc[i]==b.old_id.iloc[j]:kind='fallback';onrail=onmajor=False
 else:
  onrail=frac(shared,rt,rg)>=.3;onmajor=frac(shared,mt,mg)>=.3
  walk_cross=any(wg[k].crosses(shared) for k in wt.query(shared,predicate='intersects')) if onrail or onmajor else False
  if onrail and not(len(lt.query(shared,predicate='dwithin',distance=20)) or walk_cross):blocked['rail']+=1;continue
  inferred_edge=onmajor and not(len(pt.query(shared,predicate='dwithin',distance=20)) or walk_cross)
  if inferred_edge and frac(shared,ft,free)>=.3:blocked['major']+=1;continue
  kind='inferred_unverified' if inferred_edge else 'rail_with_crossing' if onrail else 'major_with_crossing' if onmajor else 'ordinary'
  if inferred_edge:inferred.append({'block_a':i,'block_b':j,'geometry':shared})
 add(i,j,reps[i].distance(reps[j])+(600 if kind=='inferred_unverified' else 100 if onrail or onmajor else 30),kind);newcount+=1
# Add possible but unverified crossings for primary/secondary boundaries omitted from the strict graph.
existing={(min(i,j),max(i,j)) for i,j,_ in edges};oldnodes=np.flatnonzero(active&~child);oldtree=STRtree(geo[oldnodes]);candidate_pairs=oldtree.query(geo[oldnodes],predicate='intersects').T
for ai,aj in candidate_pairs:
 i=int(oldnodes[ai]);j=int(oldnodes[aj])
 if i>=j or (i,j) in existing:continue
 shared=geo[i].boundary.intersection(geo[j].boundary)
 if shared.length<1:continue
 if frac(shared,rt,rg)>=.3 or frac(shared,ft,free)>=.3:continue
 if frac(shared,mt,mg)<.3:continue
 add(i,j,reps[i].distance(reps[j])+600,'inferred_unverified');inferred.append({'block_a':i,'block_b':j,'geometry':shared})
print('new edges',newcount,'blocked',blocked,'inferred crossings',len(inferred),flush=True)
if inferred:g.GeoDataFrame(inferred,crs=32651).to_file(out/'inferred_crossing_candidates.gpkg',driver='GPKG')
# Imagery background sampled once per node; 250m data are not treated as 250m-resolution new detail for VIIRS.
with rasterio.open('outputs/city_imagery/city_features.tif') as ds:
 names=['viirs_index','ndvi','ndbi','s1_vv_db'];bands=[ds.descriptions.index(x)+1 for x in names];xy=np.column_stack([shapely.get_x(reps),shapely.get_y(reps)]);vals=np.array(list(ds.sample(xy,indexes=bands)),dtype='float32');vals[vals==-9999]=np.nan
ntl,ndvi,ndbi,vv=vals.T;factor=np.where(np.isfinite(ntl)&(ntl<3),1.15,np.where(np.isfinite(ntl)&(ntl>21),.9,1))*np.where(np.isfinite(ndvi)&(ndvi>.35),1.8,1)*np.where(np.isfinite(vv)&(vv>-8)&np.isfinite(ndbi)&(ndbi>0),.9,1)
# Anchor mall footprints and original POI peaks. Same street block can carry multiple anchors.
mall=g.read_file('data/raw/osm/gis_osm_pois_a_free_1.shp').to_crs(32651);mall=mall[mall.fclass=='mall'];city=g.read_file('data/raw/shanghai_boundary.geojson').to_crs(32651).geometry.union_all();mall=mall[mall.geometry.representative_point().within(city)];centers=g.read_file('outputs/centers.geojson').to_crs(32651)
seed_by_node={};seed_info=[];anchor_status=[]
def locate(pt):
 ix=idx.query(pt,predicate='intersects');candidates=[int(active_ids[k]) for k in ix if geo[int(active_ids[k])].covers(pt)]
 return min(candidates,key=lambda k:geo[k].area) if candidates else None
def put(node,kind,label):
 if node not in seed_by_node:
  sid=len(seed_info)+1;seed_by_node[node]=sid;seed_info.append({'seed_id':sid,'seed_kind':kind,'anchor_labels':[],'seed_block':node})
 sid=seed_by_node[node];seed_info[sid-1]['anchor_labels'].append(label);return sid
for row in mall.itertuples():
 p=row.geometry.representative_point();node=locate(p);status='matched' if node is not None else 'unmatched_water_or_data'
 if node is not None:put(node,'commercial_aoi',str(row.name or row.osm_id))
 anchor_status.append({'osm_id':str(row.osm_id),'name':row.name,'status':status,'geometry':p})
for row in centers.itertuples():
 node=locate(row.geometry)
 if node is not None:put(node,'poi_density',str(row.id))
supplement=g.read_file(out/'additional_poi_clusters.geojson').to_crs(32651)
for row in supplement.itertuples():
 for part in shapely.get_parts(row.geometry):
  if part.geom_type!='Polygon':continue
  node=locate(part.representative_point())
  if node is not None:put(node,'poi_cluster',str(row.source_id))
print('commercial nodes',len(seed_info),'aoi matched',sum(x['status']=='matched' for x in anchor_status),'of',len(anchor_status),flush=True)
# Place local anchors on every inhabited-space bin without a commercial anchor; all are candidates, not proven inhabited.
grid={};seed_bins=set()
for node in seed_by_node:
 p=reps[node];seed_bins.add((int(np.floor(p.x/1500)),int(np.floor(p.y/1500))))
for node in np.flatnonzero(active):
 p=reps[node];key=(int(np.floor(p.x/1500)),int(np.floor(p.y/1500)))
 if key in seed_bins:continue
 cx=(key[0]+.5)*1500;cy=(key[1]+.5)*1500;d=(p.x-cx)**2+(p.y-cy)**2
 if key not in grid or d<grid[key][0]:grid[key]=(d,int(node))
for _,node in grid.values():put(node,'local','')
print('with grid local anchors',len(seed_info),flush=True)
# Every graph component must contain an anchor, even across disconnected islands / mapped barriers.
rr=[];cc=[]
for i,j,_ in edges:rr.extend((i,j));cc.extend((j,i))
graph=coo_matrix((np.ones(len(rr),dtype='int8'),(rr,cc)),shape=(n,n)).tocsr();nc,comp=connected_components(graph,directed=False)
active_nodes=np.flatnonzero(active);anchored=set(comp[list(seed_by_node)]);missing={int(comp[i]) for i in active_nodes}-anchored
area=b.geometry.area.to_numpy();order=np.argsort(comp[active_nodes],kind='stable');sorted_nodes=active_nodes[order];sorted_comp=comp[sorted_nodes];starts=np.r_[0,np.flatnonzero(np.diff(sorted_comp))+1];ends=np.r_[starts[1:],len(sorted_nodes)]
for start,end in zip(starts,ends):
 k=int(sorted_comp[start])
 if k not in missing:continue
 nodes=sorted_nodes[start:end];node=int(nodes[np.argmax(area[nodes])]);put(node,'local_component','')
print('graph components',nc,'local additional',len(missing),'total seeds',len(seed_info),flush=True)
# Multisource region growing with no arbitrary cutoff; every active component has a seed.
dist=np.full(n,np.inf);labels=np.zeros(n,dtype='int32');q=[]
for node,sid in seed_by_node.items():dist[node]=0;labels[node]=sid;heapq.heappush(q,(0.,sid,node))
while q:
 d,sid,i=heapq.heappop(q)
 if d>dist[i]+1e-7 or labels[i]!=sid:continue
 for j,cost in adj[i]:
  nd=d+cost*(factor[i]+factor[j])/2
  if nd<dist[j]:dist[j]=nd;labels[j]=sid;heapq.heappush(q,(nd,sid,j))
assert np.all(labels[active]>0)
b['unit_seed']=labels;b['cost_proxy']=np.where(np.isfinite(dist),dist,-1);b['ntl_context']=np.where(np.isfinite(ntl),ntl,-9999);b['ndvi_context']=np.where(np.isfinite(ndvi),ndvi,-9999);b.to_file(out/'activity_block_assignment.gpkg',driver='GPKG')
g.GeoDataFrame(anchor_status,crs=32651).to_crs(4326).to_file(out/'mall_anchor_audit.geojson',driver='GeoJSON');(out/'activity_seeds.json').write_text(json.dumps(seed_info,ensure_ascii=False));m={'active_land_blocks':int(active.sum()),'fallback_split_blocks':int(child.sum()),'aoi_total':len(anchor_status),'aoi_matched':sum(x['status']=='matched' for x in anchor_status),'seed_count':len(seed_info),'commercial_seed_count':sum(x['seed_kind']!='local' and x['seed_kind']!='local_component' for x in seed_info),'local_seed_count':sum(x['seed_kind'] in ('local','local_component') for x in seed_info),'unassigned_land_blocks':int((active&(labels==0)).sum()),'accepted_graph_edges':len(edges),'reused_edges':len(links),'added_split_edges':newcount,'inferred_unverified_crossing_edges':len(inferred),'blocked_new_edges':blocked,'ntl_context_blocks':int((np.isfinite(ntl)&active).sum()),'s2_context_blocks':int((np.isfinite(ndvi)&active).sum())};(out/'activity_metrics.json').write_text(json.dumps(m,ensure_ascii=False,indent=2));print(json.dumps(m,ensure_ascii=False,indent=2),flush=True)
