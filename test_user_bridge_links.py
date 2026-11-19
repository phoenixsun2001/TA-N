from pathlib import Path
import json,heapq
import geopandas as g,numpy as np,shapely
from shapely.geometry import Point
from shapely.strtree import STRtree
root=Path('outputs/user_data_test');b=g.read_file(root/'repaired_block_assignment.gpkg').to_crs(32651);geom=b.geometry.to_numpy();n=len(b);active=(~b.is_water.to_numpy())&(~b.large_flag.to_numpy());tree=STRtree(geom);reps=np.array([v.representative_point() for v in geom]);edges=g.read_file(root/'repaired_passage_links.gpkg').to_crs(32651);adj=[[] for _ in range(n)]
for e in edges.itertuples():
 i=int(e.block_a);j=int(e.block_b);adj[i].append((j,float(e.cost_m_proxy)));adj[j].append((i,float(e.cost_m_proxy)))
roads=g.read_file('data/raw/osm/gis_osm_roads_free_1.shp').to_crs(32651);roads=roads[roads.fclass.isin(['motorway','trunk','primary','secondary'])&(roads.bridge!='T')&(roads.tunnel!='T')&(roads.layer.fillna(0)==0)];rails=g.read_file('data/raw/osm/gis_osm_railways_free_1.shp').to_crs(32651);rails=rails[rails.fclass.isin(['rail','light_rail','tram'])&(rails.bridge!='T')&(rails.tunnel!='T')&(rails.layer.fillna(0)==0)];water=STRtree(geom[b.is_water.to_numpy()]);rt=STRtree(roads.geometry.to_numpy());tt=STRtree(rails.geometry.to_numpy());cand=g.read_file('data/user_shanghai/pedestrian_review_candidates.gpkg').to_crs(32651);cand=cand[cand.explicit_bridge_or_tunnel];new=[]
for ix,row in cand.iterrows():
 way=row.geometry
 if way.geom_type!='LineString' or way.length<10:continue
 if not(len(rt.query(way,predicate='intersects')) or len(tt.query(way,predicate='intersects')) or len(water.query(way,predicate='intersects'))):continue
 pts=[Point(way.coords[0]),Point(way.coords[-1])];sides=[]
 for pt in pts:
  options=[int(i) for i in tree.query(pt,predicate='dwithin',distance=35) if active[int(i)]]
  if not options:sides=[];break
  sides.append(min(options,key=lambda k:geom[k].distance(pt)))
 if len(sides)!=2 or sides[0]==sides[1]:continue
 i,j=sides
 if geom[i].boundary.intersection(geom[j].boundary).length>=1:continue
 cost=float(reps[i].distance(pts[0])+way.length+reps[j].distance(pts[1])+100);adj[i].append((j,cost));adj[j].append((i,cost));new.append({'block_a':i,'block_b':j,'kind':'user_bridge_candidate','cost_m_proxy':cost,'source_osm_id':f'user2025-row-{ix}','geometry':way})
print('eligible candidate bridge links',len(new),flush=True)
if new:g.GeoDataFrame(new,crs=32651).to_file(root/'user_bridge_candidate_links.gpkg',driver='GPKG')
seeds=json.loads((root/'road_seed_info.json').read_text());ntl=b.ntl.to_numpy();ndvi=b.ndvi.to_numpy();ndbi=b.ndbi.to_numpy();vv=b.s1_vv.to_numpy();factor=np.where(np.isfinite(ntl)&(ntl<3),1.15,np.where(np.isfinite(ntl)&(ntl>21),.9,1))*np.where(np.isfinite(ndvi)&(ndvi>.35),1.8,1)*np.where(np.isfinite(vv)&(vv>-8)&np.isfinite(ndbi)&(ndbi>0),.9,1);dist=np.full(n,np.inf);labels=np.zeros(n,dtype='int32');queue=[]
for s in seeds:
 node=int(s['seed_block']);sid=int(s['seed_id']);dist[node]=0;labels[node]=sid;heapq.heappush(queue,(0.,sid,node))
while queue:
 d,sid,i=heapq.heappop(queue)
 if d>dist[i]+1e-7 or labels[i]!=sid:continue
 for j,c in adj[i]:
  nd=d+c*(factor[i]+factor[j])/2
  if nd<dist[j]:dist[j]=nd;labels[j]=sid;heapq.heappush(queue,(nd,sid,j))
changed=active&(labels!=b.unit_seed.to_numpy());m={'explicit_bridge_tunnel_candidates':len(cand),'eligible_new_links':len(new),'changed_land_blocks':int(changed.sum()),'changed_block_area_km2':float(b.area[changed].sum()/1e6),'note':'Candidate links require current on-the-ground verification; cross-section changes are not evidence of improvement.'};(root/'bridge_candidate_test.json').write_text(json.dumps(m,indent=2));print(json.dumps(m,indent=2))
if changed.any():
 b['unit_seed']=labels;b['cost_proxy']=np.where(np.isfinite(dist),dist,-1);b.to_file(root/'repaired_block_assignment_bridge_test.gpkg',driver='GPKG')
