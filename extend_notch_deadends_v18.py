from pathlib import Path
import json,warnings
import geopandas as g,numpy as np,pandas as pd,shapely
from shapely.geometry import Point,LineString,Polygon
from shapely.ops import split
from shapely.strtree import STRtree
from convert_gcj02 import convert
out=Path('outputs/deadend_extend_v18');d=g.read_file(out/'blocks_projected.gpkg').to_crs(32651);geo=d.geometry.to_numpy();tree=STRtree(geo);roads=g.read_file('data/raw/osm/gis_osm_roads_free_1.shp').to_crs(32651);extra=g.read_file('data/user_shanghai/national_roads_candidate_filtered.gpkg').to_crs(32651);roads=g.GeoDataFrame(pd.concat([roads,extra],ignore_index=True),crs=32651);roads=roads[roads.fclass.isin(['motorway','motorway_link','trunk','trunk_link','primary','primary_link','secondary','secondary_link','tertiary','tertiary_link','residential','unclassified'])&(roads.tunnel!='T')&(roads.bridge!='T')&(roads.layer.fillna(0)==0)];candidates={}
for r in roads.itertuples():
 way=r.geometry
 if way.geom_type!='LineString' or way.length<50:continue
 for w in [way,LineString(list(way.coords)[::-1])]:
  end=Point(w.coords[-1]);prev=w.interpolate(w.length-min(20,w.length/2));v=np.array([end.x-prev.x,end.y-prev.y]);n=np.linalg.norm(v)
  if n<1e-6:continue
  v/=n;ray=LineString([(end.x,end.y),(end.x+v[0]*100,end.y+v[1]*100)])
  for j in tree.query(end,predicate='dwithin',distance=15):
   p=geo[j]
   if p.area<50000 or p.contains(end):continue
   for seg in shapely.get_parts(ray.intersection(p)):
    if seg.geom_type!='LineString' or seg.is_empty or seg.length<10:continue
    a,b=Point(seg.coords[0]),Point(seg.coords[-1]);length=end.distance(b)
    if length>100 or p.boundary.distance(a)>.01 or p.boundary.distance(b)>.01:continue
    candidates.setdefault(int(j),[]).append((length,ray,str(r.name or ''),str(r.fclass)))
rows=[];ext=[]
for i,row in d.iterrows():
 attrs=row.drop(labels='geometry').to_dict();parts=[row.geometry]
 for length,ray,name,kind in sorted(candidates.get(i,[]),key=lambda v:v[0]):
  nxt=[];changed=False
  for p in parts:
   try:q=[z for z in split(p,ray).geoms if z.geom_type=='Polygon' and z.area>=100]
   except Exception:q=[]
   if len(q)>1 and all(z.is_valid and not z.interiors for z in q) and abs(sum(z.area for z in q)-p.area)<.1:nxt.extend(q);changed=True
   else:nxt.append(p)
  if changed:
   parts=nxt;ext.append({'piece_v18':int(attrs['piece_v18']),'road_name':name,'road_class':kind,'extension_m':length,'geometry':ray.intersection(row.geometry)})
 for p in parts:rows.append({**attrs,'notch_deadend_split':len(parts)>1,'geometry':p})
res=g.GeoDataFrame(rows,crs=32651);res['piece_final']=np.arange(1,len(res)+1);res['area_km2']=res.area/1e6;res.to_file(out/'blocks_projected_final.gpkg',driver='GPKG')
if ext:g.GeoDataFrame(ext,crs=32651).to_file(out/'notch_extensions.gpkg',driver='GPKG')
for frame,name in [(res,'blocks_gcj02_final.geojson')]+([(g.GeoDataFrame(ext,crs=32651),'notch_extensions_gcj02.geojson')] if ext else []):
 x=frame.to_crs(4326);x.geometry=shapely.transform(x.geometry.array,convert,include_z=False);bad=~x.is_valid;x.loc[bad,'geometry']=shapely.make_valid(x.loc[bad,'geometry'].array);x=x.explode(index_parts=False,ignore_index=True)
 if name.startswith('blocks'):x=x[x.geom_type=='Polygon'];x.geometry=[Polygon(z.exterior) for z in x.geometry]
 x.crs=None
 with warnings.catch_warnings():warnings.filterwarnings('ignore',message="'crs' was not provided.*");x.to_file(out/name,driver='GeoJSON')
 if name.startswith('blocks'):gcj={'features':len(x),'invalid':int((~x.is_valid).sum()),'holes':sum(len(z.interiors) for z in x.geometry)}
m={'source_blocks':len(d),'candidates':sum(map(len,candidates.values())),'notch_extensions':len(ext),'final_blocks':len(res),'area_delta_m2':float(res.area.sum()-d.area.sum()),'invalid':int((~res.is_valid).sum()),'holes':sum(len(z.interiors) for z in res.geometry),'gcj':gcj};(out/'notch_metrics.json').write_text(json.dumps(m,indent=2));print(json.dumps(m,indent=2))
