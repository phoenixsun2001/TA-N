from pathlib import Path
import json,warnings
import geopandas as g,numpy as np,pandas as pd,shapely
from shapely.geometry import Point,LineString
from shapely.ops import split,linemerge,nearest_points
from shapely.strtree import STRtree
from convert_gcj02 import convert
out=Path('outputs/deadend_extend_v18');out.mkdir(exist_ok=True)
d=g.read_file('outputs/poi_direct_v17/blocks_projected.gpkg').to_crs(32651);r=g.read_file('data/raw/osm/gis_osm_roads_free_1.shp').to_crs(32651);extra=g.read_file('data/user_shanghai/national_roads_candidate_filtered.gpkg').to_crs(32651);r=g.GeoDataFrame(pd.concat([r,extra],ignore_index=True),crs=32651)
classes=['motorway','motorway_link','trunk','trunk_link','primary','primary_link','secondary','secondary_link','tertiary','tertiary_link','residential','unclassified'];r=r[r.fclass.isin(classes)&(r.tunnel!='T')&(r.bridge!='T')&(r.layer.fillna(0)==0)].reset_index(drop=True);rg=r.geometry.to_numpy();rt=STRtree(rg)
rows=[];cuts=[];extended=[];stats={'candidate_ways':0,'too_far_or_no_forward_boundary':0,'no_valid_split':0};
for idx,row in d.iterrows():
 attrs=row.drop(labels='geometry').to_dict();area=row.geometry;parts=[area]
 if area.area>5000:
  near=rt.query(area,predicate='intersects');clipped=[z for j in near for z in shapely.get_parts(rg[j].intersection(area)) if z.geom_type=='LineString' and z.length>1]
  try:network=linemerge(shapely.union_all(clipped)) if len(clipped)>1 else clipped[0] if clipped else None
  except Exception:network=None
  if network is not None:
   ways=sorted([z for z in shapely.get_parts(network) if z.geom_type=='LineString' and z.length>=50],key=lambda z:-z.length)
   for way in ways:
    a,b=Point(way.coords[0]),Point(way.coords[-1]);da,db=area.boundary.distance(a),area.boundary.distance(b)
    if da<=2 and db>2:forward=way
    elif db<=2 and da>2:forward=LineString(list(way.coords)[::-1]);a,b=b,a;da,db=db,da
    else:continue
    stats['candidate_ways']+=1;end=b;prev=forward.interpolate(max(0,forward.length-min(20,forward.length/2)));vec=np.array([end.x-prev.x,end.y-prev.y]);norm=np.linalg.norm(vec)
    if norm<1e-6:continue
    vec/=norm;ray=LineString([(end.x,end.y),(end.x+vec[0]*100,end.y+vec[1]*100)]);hit=ray.intersection(area.boundary);points=[]
    for obj in shapely.get_parts(hit):
     if obj.geom_type=='Point' and obj.distance(end)>.1:points.append(obj)
     elif obj.geom_type=='MultiPoint':points.extend([v for v in obj.geoms if v.distance(end)>.1])
    if not points:stats['too_far_or_no_forward_boundary']+=1;continue
    tip=min(points,key=lambda z:z.distance(end));delta=tip.distance(end)
    if delta>100:stats['too_far_or_no_forward_boundary']+=1;continue
    start=nearest_points(a,area.boundary)[1];startvec=np.array([start.x-a.x,start.y-a.y]);vnorm=np.linalg.norm(startvec)
    if vnorm<1e-6:startvec=np.array([a.x-forward.coords[1][0],a.y-forward.coords[1][1]]);vnorm=np.linalg.norm(startvec)
    if vnorm<1e-6:continue
    begin=np.array([start.x,start.y])+startvec/vnorm*.05;finish=np.array([tip.x,tip.y])+vec*.05;full=LineString([tuple(begin),*forward.coords,tuple(finish)])
    new=[];changed=False
    for part in parts:
     try:q=[v for v in split(part,full).geoms if v.geom_type=='Polygon' and v.area>=100]
     except Exception:q=[]
     if len(q)>1 and all(v.is_valid and not v.interiors for v in q) and abs(sum(v.area for v in q)-part.area)<.1:new.extend(q);changed=True
     else:new.append(part)
    if changed:
     parts=new;mid=forward.interpolate(.5,normalized=True);nearest=int(near[np.argmin(shapely.distance(rg[near],mid))]);rec={'direct_piece_id':int(attrs['direct_piece_id']),'road_name':str(r.name.iloc[nearest] or ''),'road_class':str(r.fclass.iloc[nearest]),'extension_m':float(delta),'geometry':LineString([end,tip])};extended.append(rec);cuts.append({**rec,'geometry':full.intersection(area)})
    else:stats['no_valid_split']+=1
 for p in parts:rows.append({**attrs,'deadend_split':len(parts)>1,'geometry':p})
res=g.GeoDataFrame(rows,crs=32651);res['piece_v18']=np.arange(1,len(res)+1);res['area_km2']=res.area/1e6;res.to_file(out/'blocks_projected.gpkg',driver='GPKG')
for dat,name in [(extended,'extensions'),(cuts,'cutlines')]:
 if dat:g.GeoDataFrame(dat,crs=32651).to_file(out/(name+'.gpkg'),driver='GPKG')
for frame,name in [(res,'blocks_gcj02.geojson')]+([(g.GeoDataFrame(extended,crs=32651),'extensions_gcj02.geojson')] if extended else [])+([(g.GeoDataFrame(cuts,crs=32651),'cutlines_gcj02.geojson')] if cuts else []):
 x=frame.to_crs(4326);x.geometry=shapely.transform(x.geometry.array,convert,include_z=False);bad=~x.is_valid;x.loc[bad,'geometry']=shapely.make_valid(x.loc[bad,'geometry'].array);x=x.explode(index_parts=False,ignore_index=True)
 if name.startswith('blocks'):
  x=x[x.geom_type=='Polygon'];from shapely.geometry import Polygon
  x.geometry=[Polygon(z.exterior) for z in x.geometry]
 x.crs=None
 with warnings.catch_warnings():warnings.filterwarnings('ignore',message="'crs' was not provided.*");x.to_file(out/name,driver='GeoJSON')
 if name.startswith('blocks'):gcj={'features':len(x),'invalid':int((~x.is_valid).sum()),'holes':sum(len(z.interiors) for z in x.geometry)}
m={'input_blocks':len(d),'output_blocks':len(res),'source_blocks_split':int(sum(res.groupby('direct_piece_id').size()>1)),'extensions':len(extended),'extension_max_m':max([x['extension_m'] for x in extended],default=0),'extension_median_m':float(np.median([x['extension_m'] for x in extended])) if extended else None,'stats':stats,'area_change_m2':float(res.area.sum()-d.area.sum()),'invalid':int((~res.is_valid).sum()),'holes':sum(len(x.interiors) for x in res.geometry),'gcj':gcj};(out/'metrics.json').write_text(json.dumps(m,indent=2));print(json.dumps(m,indent=2),flush=True)
