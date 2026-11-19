from pathlib import Path
import json,geopandas as g,shapely
from shapely.geometry import Polygon
from shapely.ops import triangulate
out=Path('outputs/revised_poc');parent=g.read_file(out/'activity_units_parent.gpkg').to_crs(32651)
def polys(x):
 if x.geom_type=='Polygon':return [x]
 if hasattr(x,'geoms'):return [p for q in x.geoms for p in polys(q)]
 return []
def noholes(p):
 if not p.interiors:return [p]
 x0,y0,x1,y1=p.bounds;cuts=sorted(set(Polygon(r).representative_point().y for r in p.interiors));bounds=[y0-1]+cuts+[y1+1];result=[]
 for low,high in zip(bounds[:-1],bounds[1:]):
  for q in polys(p.intersection(shapely.box(x0-1,low,x1+1,high))):
   if q.area<=1e-9:continue
   if q.interiors:
    for t in triangulate(q):
     for x in polys(t.intersection(q)):
      if x.area>1e-9:
       if x.interiors:raise RuntimeError('ring remains')
       result.append(x)
   else:result.append(q)
 return result
rows=[];rings=0
for row in parent.itertuples():
 old=polys(row.geometry);rings+=sum(len(v.interiors) for v in old);pieces=[]
 for v in old:pieces.extend(noholes(v))
 for i,v in enumerate(pieces):rows.append({'component_id':f'{row.unit_seed}-{i+1:04d}','unit_seed':int(row.unit_seed),'zone_type':row.zone_type,'fallback_split':row.fallback_split,'inferred_crossing':row.inferred_crossing,'hole_cut':len(pieces)>len(old),'geometry':v})
f=g.GeoDataFrame(rows,crs=32651);print('pieces',len(f),'invalid',int((~f.is_valid).sum()),'holes',sum(len(x.interiors) for x in f.geometry),flush=True);assert f.is_valid.all() and all(len(x.interiors)==0 for x in f.geometry)
f.to_file(out/'activity_components_no_holes.gpkg',driver='GPKG');(out/'components_metrics.json').write_text(json.dumps({'parent_units':len(parent),'components':len(f),'original_interior_rings':rings,'invalid':int((~f.is_valid).sum()),'remaining_interior_rings':sum(len(x.interiors) for x in f.geometry)},indent=2));print('done',flush=True)
