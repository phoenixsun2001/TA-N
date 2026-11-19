from pathlib import Path
import json,warnings
import geopandas as g
import numpy as np,shapely
from shapely.geometry import Polygon,LineString
from shapely.ops import nearest_points,split
from convert_gcj02 import convert
out=Path('outputs/closed_small_holes_v8');out.mkdir(parents=True,exist_ok=True)
d=g.read_file('outputs/improve_v7/blocks_projected.gpkg').to_crs(32651)
rows=[];cuts=[];failed=[];small_before=0
for row in d.itertuples(index=False):
 attrs=row._asdict();p=attrs.pop('geometry');work=[p];finished=[];attempted=set()
 while work:
  p=work.pop();rings=[r for r in p.interiors if Polygon(r).area<=1000]
  if not rings:finished.append(p);continue
  ring=min(rings,key=lambda x:Polygon(x).area);a,b=nearest_points(ring,p.exterior);small_before+=1
  key=(round(a.x,4),round(a.y,4),round(b.x,4),round(b.y,4))
  length=a.distance(b)
  if length<1e-8 or key in attempted:
   finished.append(p);failed.append(int(attrs['block_id']));continue
  attempted.add(key);dx=(b.x-a.x)/length;dy=(b.y-a.y)/length;eps=.01
  line=LineString([(a.x-dx*eps,a.y-dy*eps),(b.x+dx*eps,b.y+dy*eps)])
  try:parts=[z for z in split(p,line).geoms if z.geom_type=='Polygon' and z.area>0]
  except Exception:parts=[]
  # A radial cut can open a hole without creating two polygons; polygon splitting
  # needs two cuts, so use a second cut at the opposite end of the same ring.
  if len(parts)<2:
   coords=np.asarray(ring.coords);dist=np.hypot(coords[:,0]-a.x,coords[:,1]-a.y);far=shapely.Point(coords[np.argmax(dist)]);c,e=nearest_points(far,p.exterior);ll=c.distance(e)
   if ll>1e-8:
    dx2=(e.x-c.x)/ll;dy2=(e.y-c.y)/ll;line2=LineString([(c.x-dx2*eps,c.y-dy2*eps),(e.x+dx2*eps,e.y+dy2*eps)])
    try:
     from shapely.ops import polygonize
     network=shapely.union_all([p.boundary,line.intersection(p),line2.intersection(p)])
     parts=[z for z in polygonize(network) if p.covers(z.representative_point()) and z.area>0]
    except Exception:parts=[]
    if len(parts)>=2:cuts.append({'block_id':int(attrs['block_id']),'geometry':line2.intersection(p)})
  if len(parts)<2 or abs(sum(z.area for z in parts)-p.area)>.05:
   center=Polygon(ring).representative_point();near=nearest_points(center,p.exterior)[1];span=max(p.bounds[2]-p.bounds[0],p.bounds[3]-p.bounds[1])*3+1
   distance=center.distance(near)
   if distance>1e-8:
    ux=(near.x-center.x)/distance;uy=(near.y-center.y)/distance
    full=LineString([(center.x-ux*span,center.y-uy*span),(center.x+ux*span,center.y+uy*span)])
    try:parts=[z for z in split(p,full).geoms if z.geom_type=='Polygon' and z.area>0]
    except Exception:parts=[]
    if len(parts)>=2 and abs(sum(z.area for z in parts)-p.area)<=.05:line=full
   if len(parts)<2 or abs(sum(z.area for z in parts)-p.area)>.05:
    finished.append(p);failed.append(int(attrs['block_id']));continue
  cuts.append({'block_id':int(attrs['block_id']),'geometry':line.intersection(p)})
  work.extend(parts)
 for n,z in enumerate(finished):rows.append({**attrs,'part':n+1,'geometry':z})
res=g.GeoDataFrame(rows,crs=32651);res.insert(0,'piece_id',np.arange(1,len(res)+1));res['area_km2']=res.area/1e6;res['holes_remaining']=[len(p.interiors) for p in res.geometry];res.to_file(out/'blocks_projected.gpkg',driver='GPKG')
if cuts:g.GeoDataFrame(cuts,crs=32651).to_file(out/'added_cutlines.gpkg',driver='GPKG')
x=res.to_crs(4326);x.geometry=shapely.transform(x.geometry.array,convert,include_z=False);bad=~x.is_valid;x.loc[bad,'geometry']=shapely.make_valid(x.loc[bad,'geometry'].array);x=x.explode(index_parts=False,ignore_index=True);x=x[x.geom_type=='Polygon'];x.crs=None
with warnings.catch_warnings():
 warnings.filterwarnings('ignore',message="'crs' was not provided.*");x.to_file(out/'blocks_gcj02.geojson',driver='GeoJSON')
remaining=sum(Polygon(r).area<=1000 for p in res.geometry for r in p.interiors)
m={'original_blocks':len(d),'output_parts':len(res),'small_holes_before':sum(Polygon(r).area<=1000 for p in d.geometry for r in p.interiors),'small_holes_remaining':remaining,'failed_source_blocks':len(set(failed)),'cutlines':len(cuts),'projected_invalid':int((~res.is_valid).sum()),'area_change_m2':float(res.area.sum()-d.area.sum()),'gcj_features':len(x),'gcj_invalid':int((~x.is_valid).sum())};(out/'metrics.json').write_text(json.dumps(m,indent=2));print(json.dumps(m,indent=2))
