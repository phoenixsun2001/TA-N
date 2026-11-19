from pathlib import Path
import json
import geopandas as g, shapely
from shapely.geometry import Polygon
from shapely.ops import triangulate
out=Path('outputs/revised_poc/gcj02');path=out/'activity_components_no_holes_gcj02.geojson';df=g.read_file(path)
def parts(x):
 if x.geom_type=='Polygon':return [x]
 if hasattr(x,'geoms'):return [p for q in x.geoms for p in parts(q)]
 return []
def noholes(p):
 if not p.interiors:return [p]
 x0,y0,x1,y1=p.bounds;mx=max((x1-x0)*.01,1e-5);cuts=sorted(set(Polygon(r).representative_point().y for r in p.interiors));ys=[y0-mx]+cuts+[y1+mx];result=[]
 for a,b in zip(ys[:-1],ys[1:]):
  for q in parts(p.intersection(shapely.box(x0-mx,a,x1+mx,b))):
   if q.area<=1e-24:continue
   if q.interiors:
    for t in triangulate(q):
     for v in parts(t.intersection(q)):
      if v.area>1e-24:
       if v.interiors:raise RuntimeError('unresolved internal ring')
       result.append(v)
   else:result.append(q)
 return result
rows=[];discarded=[]
for row in df.itertuples(index=False):
 attrs=row._asdict();geometry=attrs.pop('geometry');polys=parts(geometry)
 if not polys:discarded.append(str(attrs.get('component_id')));continue
 ps=[]
 for p in polys:ps.extend(noholes(p))
 for i,p in enumerate(ps):
  if not p.is_valid:
   ps2=parts(shapely.make_valid(p))
   for k,v in enumerate(ps2):
    if v.interiors:raise RuntimeError('invalid repair created internal ring')
    rows.append({**attrs,'gcj_part':f'{i+1}.{k+1}','geometry':v})
  else:rows.append({**attrs,'gcj_part':str(i+1),'geometry':p})
f=g.GeoDataFrame(rows,crs=None);assert f.is_valid.all() and (f.geom_type=='Polygon').all() and sum(len(x.interiors) for x in f.geometry)==0
f.to_file(path,driver='GeoJSON');(out/'components_repair_metrics.json').write_text(json.dumps({'input_features':len(df),'output_polygon_features':len(f),'collapsed_to_line_count':len(discarded),'collapsed_component_ids':discarded,'remaining_invalid':int((~f.is_valid).sum()),'remaining_holes':sum(len(x.interiors) for x in f.geometry)},ensure_ascii=False,indent=2));print('input',len(df),'output',len(f),'collapsed',len(discarded),flush=True)
