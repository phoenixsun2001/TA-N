from pathlib import Path
import json,warnings
import geopandas as g,numpy as np,shapely
from shapely.geometry import Point,LineString
from shapely.ops import split
from shapely.strtree import STRtree
from convert_gcj02 import convert
out=Path('outputs/riverfront_v15');out.mkdir(exist_ok=True)
d=g.read_file('outputs/deadend_notches_v14/blocks_projected.gpkg').to_crs(32651)
water=g.read_file('data/raw/osm/gis_osm_water_a_free_1.shp').to_crs(32651);water=water[(water.fclass=='riverbank')|((water.fclass=='water')&(water.area>=100000))];waters=shapely.make_valid(water.geometry.array);wt=STRtree(waters)
roads=g.read_file('data/raw/osm/gis_osm_roads_free_1.shp').to_crs(32651);extra=g.read_file('data/user_shanghai/national_roads_candidate_filtered.gpkg').to_crs(32651)
import pandas as pd
roads=g.GeoDataFrame(pd.concat([roads,extra],ignore_index=True),crs=32651);types=['motorway','motorway_link','trunk','trunk_link','primary','primary_link','secondary','secondary_link','tertiary','tertiary_link','residential','unclassified'];roads=roads[roads.fclass.isin(types)&(roads.tunnel!='T')&(roads.layer.fillna(0)==0)&(roads.bridge!='T')].copy().reset_index(drop=True);rg=roads.geometry.to_numpy();rt=STRtree(rg)
rows=[];lines=[];review=[];count=0
for row in d.itertuples(index=False):
 attrs=row._asdict();geom=attrs.pop('geometry');box=geom.minimum_rotated_rectangle;xy=np.asarray(box.exterior.coords);edges=np.linalg.norm(np.diff(xy,axis=0),axis=1);long=max(edges);short=min(edges);nearwater=wt.query(geom.boundary,predicate='dwithin',distance=50)
 if geom.area<100000 or short<1 or long<800 or long/short<4 or len(nearwater)==0:
  rows.append({**attrs,'river_review':False,'river_split':False,'geometry':geom});continue
 count+=1;parts=[geom];selected=[]
 candidates=rt.query(geom,predicate='intersects')
 if len(candidates):
  from shapely.ops import nearest_points,linemerge
  clipped=[part for k in candidates for part in shapely.get_parts(rg[k].intersection(geom)) if part.geom_type=='LineString' and part.length>1]
  try:merged=linemerge(shapely.union_all(clipped)) if len(clipped)>1 else clipped[0] if clipped else None
  except Exception:merged=None
  if merged is not None:
   for piece in shapely.get_parts(merged):
    if piece.geom_type!='LineString' or piece.length<30:continue
    a=Point(piece.coords[0]);b=Point(piece.coords[-1]);shore=min(min(waters[t].distance(a),waters[t].distance(b)) for t in nearwater)
    da=geom.boundary.distance(a);db=geom.boundary.distance(b)
    if shore>50 or da>20 or db>20:continue
    aa=nearest_points(a,geom.boundary)[1];bb=nearest_points(b,geom.boundary)[1]
    coords=[(aa.x,aa.y),*piece.coords,(bb.x,bb.y)]
    # Extend the short endpoint connectors through the boundary for robust splitting.
    first=Point(coords[1]);last=Point(coords[-2]);va=np.array(coords[0])-np.array([first.x,first.y]);vb=np.array(coords[-1])-np.array([last.x,last.y]);na=np.linalg.norm(va);nb=np.linalg.norm(vb)
    if na>1e-6:coords[0]=tuple(np.array(coords[0])+va/na*.05)
    if nb>1e-6:coords[-1]=tuple(np.array(coords[-1])+vb/nb*.05)
    cut=LineString(coords);nextparts=[];changed=False
    for p in parts:
     try:q=[v for v in split(p,cut).geoms if v.geom_type=='Polygon' and v.area>1]
     except Exception:q=[]
     if len(q)>=2 and min(v.area for v in q)>=5000 and abs(sum(v.area for v in q)-p.area)<.1:nextparts.extend(q);changed=True
     else:nextparts.append(p)
    if changed:
     parts=nextparts;k=int(candidates[np.argmin(shapely.distance(rg[candidates],piece.interpolate(.5,normalized=True)))]);selected.append((k,cut))
 for k,piece in selected:lines.append({'block_id':int(attrs['block_id']),'road_class':str(roads.fclass.iloc[k]),'road_name':str(roads.name.iloc[k] or ''),'geometry':piece})
 review.append({'block_id':int(attrs['block_id']),'length_m':long,'ratio':long/short,'split_count':len(parts),'geometry':geom})
 for p in parts:rows.append({**attrs,'river_review':True,'river_split':len(parts)>1,'geometry':p})
r=g.GeoDataFrame(rows,crs=32651);r['river_piece_id']=np.arange(1,len(r)+1);r['area_km2']=r.area/1e6;r.to_file(out/'blocks_projected.gpkg',driver='GPKG');g.GeoDataFrame(review,crs=32651).to_file(out/'riverfront_review.gpkg',driver='GPKG')
if lines:g.GeoDataFrame(lines,crs=32651).to_file(out/'river_split_lines.gpkg',driver='GPKG')
x=r.to_crs(4326);x.geometry=shapely.transform(x.geometry.array,convert,include_z=False);bad=~x.is_valid;x.loc[bad,'geometry']=shapely.make_valid(x.loc[bad,'geometry'].array);x=x.explode(index_parts=False,ignore_index=True);x=x[x.geom_type=='Polygon'];x.crs=None
with warnings.catch_warnings():
 warnings.filterwarnings('ignore',message="'crs' was not provided.*");x.to_file(out/'blocks_gcj02.geojson',driver='GeoJSON')
m={'input_features':len(d),'review_candidates':count,'split_source_blocks':sum(v['split_count']>1 for v in review),'split_lines':len(lines),'output_features':len(r),'area_delta_m2':float(r.area.sum()-d.area.sum()),'invalid':int((~r.is_valid).sum()),'holes':sum(len(p.interiors) for p in r.geometry),'gcj_features':len(x),'gcj_invalid':int((~x.is_valid).sum()),'gcj_holes':sum(len(p.interiors) for p in x.geometry)};(out/'metrics.json').write_text(json.dumps(m,indent=2));print(json.dumps(m,indent=2),flush=True)
