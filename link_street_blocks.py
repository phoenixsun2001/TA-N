from pathlib import Path
import json,heapq
import numpy as np,geopandas as gpd,shapely
from shapely.geometry import Point,Polygon
from shapely.ops import unary_union
from shapely.strtree import STRtree
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
out=Path('outputs/block_fit')
def read(name):return gpd.read_file(f'data/raw/osm/gis_osm_{name}_free_1.shp').to_crs(32651)
roads=read('roads');rail=read('railways');traffic=read('traffic');centers=gpd.read_file('outputs/centers.geojson').to_crs(32651)
major=roads[roads.fclass.isin(['motorway','trunk','primary','secondary'])&(roads.bridge!='T')&(roads.tunnel!='T')&(roads.layer.fillna(0)==0)]
rail=rail[rail.fclass.isin(['rail','light_rail','tram'])&(rail.bridge!='T')&(rail.tunnel!='T')&(rail.layer.fillna(0)==0)]
walk=roads[roads.fclass.isin(['footway','path','pedestrian','steps'])];bridge=walk[(walk.bridge=='T')|(walk.tunnel=='T')]
ped=traffic[traffic.fclass=='pedestrian_crossing'];level=traffic[traffic.fclass=='railway_crossing']
for name in ['xujiahui','lujiazui','jinqiao']:
 blocks=gpd.read_file(out/f'{name}_blocks.gpkg').to_crs(32651);geom=blocks.geometry.to_numpy();m=json.loads((out/f'{name}_input_metrics.json').read_text());x0,y0,x1,y1=m['domain_bounds']
 def sub(df):return df.cx[x0:x1,y0:y1]
 major2=sub(major);rail2=sub(rail);ped2=sub(ped);level2=sub(level);walk2=sub(walk);bridge2=sub(bridge);seeds=sub(centers)
 majorbuf=major2.geometry.union_all().buffer(3);railbuf=rail2.geometry.union_all().buffer(3) if len(rail2) else None
 ptree=STRtree(ped2.geometry.to_numpy());ltree=STRtree(level2.geometry.to_numpy());wtree=STRtree(walk2.geometry.to_numpy());btree=STRtree(geom)
 reps=np.array([p.representative_point() for p in geom]); pairs=btree.query(geom,predicate='intersects').T
 adj=[[] for _ in geom]; records=[];count={'ordinary':0,'major_with_crossing':0,'rail_with_crossing':0,'blocked_major':0,'blocked_rail':0,'bridge_link':0}
 def linked(i,j,cost,kind):
  adj[i].append((j,cost));adj[j].append((i,cost));records.append({'a':i,'b':j,'kind':kind,'cost_m_proxy':cost,'geometry':shapely.shortest_line(reps[i],reps[j])})
 for i,j in pairs:
  i=int(i);j=int(j)
  if i>=j or blocks.is_water[i] or blocks.is_water[j]:continue
  shared=geom[i].boundary.intersection(geom[j].boundary)
  if shared.length<1:continue
  onrail=railbuf is not None and shared.intersection(railbuf).length>=min(5,shared.length*.3)
  onmajor=shared.intersection(majorbuf).length>=min(5,shared.length*.3)
  crosslevel=len(ltree.query(shared,predicate='dwithin',distance=20))>0 if onrail else False
  crossped=len(ptree.query(shared,predicate='dwithin',distance=20))>0 if onmajor else False
  # A mapped pedestrian way that geometrically crosses the common boundary is a crossing proxy.
  crossing_walk=any(walk2.geometry.iloc[k].crosses(shared) for k in wtree.query(shared,predicate='intersects'))
  if onrail and not (crosslevel or crossing_walk):count['blocked_rail']+=1;continue
  if onmajor and not (crossped or crossing_walk):count['blocked_major']+=1;continue
  kind='rail_with_crossing' if onrail else 'major_with_crossing' if onmajor else 'ordinary';count[kind]+=1
  linked(i,j,float(reps[i].distance(reps[j]))+(100 if onmajor or onrail else 30),kind)
 # Explicit mapped pedestrian bridge/tunnel: connect land faces at opposite endpoints only if separated by >= 10 m.
 for way in bridge2.itertuples():
  if way.geometry.geom_type!='LineString' or way.geometry.length<10:continue
  endpoints=[Point(way.geometry.coords[0]),Point(way.geometry.coords[-1])];sides=[]
  for point in endpoints:
   candidates=[int(i) for i in btree.query(point,predicate='dwithin',distance=35) if not blocks.is_water[int(i)]]
   if not candidates:sides=[];break
   sides.append(min(candidates,key=lambda i:geom[i].distance(point)))
  if len(sides)==2 and sides[0]!=sides[1]:
   i,j=sides
   if geom[i].boundary.intersection(geom[j].boundary).length>=1:continue
   linked(i,j,float(reps[i].distance(endpoints[0])+way.geometry.length+reps[j].distance(endpoints[1])+100),'bridge_link');count['bridge_link']+=1
 # Map global seeds to containing / nearest land block.
 seednodes=[]
 for r in seeds.itertuples():
  candidates=[int(i) for i in btree.query(r.geometry,predicate='intersects') if not blocks.is_water[int(i)]]
  if not candidates:
   candidates=[int(i) for i in btree.query(r.geometry,predicate='dwithin',distance=100) if not blocks.is_water[int(i)]]
  if candidates:seednodes.append((int(r.id),min(candidates,key=lambda i:geom[i].distance(r.geometry))))
 dist=np.full(len(blocks),np.inf);labels=np.zeros(len(blocks),int);heap=[]
 for id,node in seednodes:
  if dist[node]==0:continue
  dist[node]=0;labels[node]=id;heapq.heappush(heap,(0,id,node))
 while heap:
  d,id,n=heapq.heappop(heap)
  if d>dist[n]+.001 or labels[n]!=id:continue
  for j,w in adj[n]:
   dd=d+w
   if dd<=2000 and dd<dist[j]:dist[j]=dd;labels[j]=id;heapq.heappush(heap,(dd,id,j))
 blocks['seed_id']=labels;blocks['reachable']=labels>0;blocks['cost_m_proxy']=np.where(np.isfinite(dist),dist,-1);blocks.to_file(out/f'{name}_block_assignment.gpkg',driver='GPKG')
 if records:gpd.GeoDataFrame(records,crs=32651).to_file(out/f'{name}_passage_links.gpkg',driver='GPKG')
 # Report resulting components and enclosed holes before topology normalization.
 comp=holes=0
 for id in sorted(set(labels)-{0}):
  merged=unary_union([geom[i] for i in range(len(geom)) if labels[i]==id]);parts=list(merged.geoms) if merged.geom_type=='MultiPolygon' else [merged];comp+=len(parts);holes+=sum(len(p.interiors) for p in parts if p.geom_type=='Polygon')
 result={'site':name,'blocks':len(blocks),'land_blocks':int((~blocks.is_water).sum()),'seeds_matched':len(seednodes),'assigned_blocks':int((labels>0).sum()),'unassigned_land_blocks':int(((labels==0)&(~blocks.is_water)).sum()),'graph_edges':count,'seed_polygon_components':comp,'interior_rings_before_normalization':holes}
 (out/f'{name}_graph_metrics.json').write_text(json.dumps(result,indent=2));print(json.dumps(result),flush=True)
