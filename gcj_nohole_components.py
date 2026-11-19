from pathlib import Path
import geopandas as g,shapely,json,warnings
from shapely.geometry import Polygon
from shapely.ops import triangulate
out=Path('outputs/revised_poc/gcj02');parent=g.read_file(out/'activity_units_parent_gcj02.geojson')
def polys(x):
 if x.geom_type=='Polygon':return [x]
 if hasattr(x,'geoms'):return [p for q in x.geoms for p in polys(q)]
 return []
def noholes(p):
 if not p.interiors:return [p]
 x0,y0,x1,y1=p.bounds;margin=max((x1-x0)*.01,1e-5);cuts=sorted(set(Polygon(r).representative_point().y for r in p.interiors));bounds=[y0-margin]+cuts+[y1+margin];res=[]
 for low,high in zip(bounds[:-1],bounds[1:]):
  for q in polys(p.intersection(shapely.box(x0-margin,low,x1+margin,high))):
   if q.area<=1e-24:continue
   if q.interiors:
    for t in triangulate(q):
     for v in polys(t.intersection(q)):
      if v.area>1e-24:
       if v.interiors:raise RuntimeError('ring remains')
       res.append(v)
   else:res.append(q)
 return res
rows=[];ring_count=0
for row in parent.itertuples(index=False):
 old=polys(row.geometry);ring_count+=sum(len(v.interiors) for v in old);pieces=[]
 for p in old:pieces.extend(noholes(p))
 for i,p in enumerate(pieces,1):
  rows.append({'component_id':f'{row.unit_seed}-{i:04d}','unit_seed':int(row.unit_seed),'zone_type':row.zone_type,'fallback_split':row.fallback_split,'inferred_crossing':row.inferred_crossing,'hole_cut':len(pieces)>len(old),'geometry':p})
f=g.GeoDataFrame(rows,crs=None);bad=~f.is_valid
if bad.any():raise RuntimeError(f'{bad.sum()} invalid components; inspect before export')
assert (f.geom_type=='Polygon').all() and all(not len(v.interiors) for v in f.geometry)
with warnings.catch_warnings():
 warnings.filterwarnings('ignore',message="'crs' was not provided.*");f.to_file(out/'activity_components_no_holes_gcj02.geojson',driver='GeoJSON')
m={'parent_units':len(parent),'components':len(f),'source_holes_cut':ring_count,'invalid':int(bad.sum()),'remaining_holes':sum(len(x.interiors) for x in f.geometry),'note':'Components are pieces of grouped units and cutting can introduce artificial boundaries.'};(out/'components_topology_metrics.json').write_text(json.dumps(m,indent=2));print(json.dumps(m,indent=2),flush=True)
