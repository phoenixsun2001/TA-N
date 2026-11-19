from pathlib import Path
import json,warnings
import geopandas as g,numpy as np,shapely
from shapely.geometry import Polygon,LineString
from shapely.ops import split
from convert_gcj02 import convert
out=Path('outputs/closed_small_holes_v8');d=g.read_file(out/'blocks_projected.gpkg');rows=[];cuts=[];failed=[]
for row in d.itertuples(index=False):
 attrs=row._asdict();p=attrs.pop('geometry');todo=[p]
 while todo:
  p=todo.pop();holes=[r for r in p.interiors if Polygon(r).area<=1000]
  if not holes:rows.append({**attrs,'geometry':p});continue
  hole=Polygon(holes[0]);center=hole.representative_point();span=max(p.bounds[2]-p.bounds[0],p.bounds[3]-p.bounds[1])*3+1;parts=[];line=None
  for angle in np.linspace(0,np.pi,24,endpoint=False):
   ux=np.cos(angle);uy=np.sin(angle);candidate=LineString([(center.x-ux*span,center.y-uy*span),(center.x+ux*span,center.y+uy*span)])
   try:q=[z for z in split(p,candidate).geoms if z.geom_type=='Polygon' and z.area>1e-10]
   except Exception:continue
   new=sum(Polygon(r).area<=1000 for z in q for r in z.interiors)
   if len(q)>1 and new<len(holes) and abs(sum(z.area for z in q)-p.area)<.05:parts=q;line=candidate;break
  if not parts:rows.append({**attrs,'geometry':p});failed.append(int(attrs['block_id']));continue
  cuts.append({'block_id':int(attrs['block_id']),'geometry':line.intersection(p)});todo.extend(parts)
r=g.GeoDataFrame(rows,crs=32651);r['piece_id']=np.arange(1,len(r)+1);r['area_km2']=r.area/1e6;r['holes_remaining']=[len(z.interiors) for z in r.geometry];r.to_file(out/'blocks_projected_final.gpkg',driver='GPKG')
if cuts:g.GeoDataFrame(cuts,crs=32651).to_file(out/'additional_cutlines.gpkg',driver='GPKG')
x=r.to_crs(4326);x.geometry=shapely.transform(x.geometry.array,convert,include_z=False);bad=~x.is_valid;x.loc[bad,'geometry']=shapely.make_valid(x.loc[bad,'geometry'].array);x=x.explode(index_parts=False,ignore_index=True);x=x[x.geom_type=='Polygon'];x.crs=None
with warnings.catch_warnings():
 warnings.filterwarnings('ignore',message="'crs' was not provided.*");x.to_file(out/'blocks_gcj02_final.geojson',driver='GeoJSON')
m={'features':len(r),'small_holes_remaining':sum(Polygon(z).area<=1000 for p in r.geometry for z in p.interiors),'failed_blocks':len(set(failed)),'additional_cutlines':len(cuts),'invalid':int((~r.is_valid).sum()),'area_change_m2':float(r.area.sum()-d.area.sum()),'gcj_features':len(x),'gcj_invalid':int((~x.is_valid).sum())};(out/'final_metrics.json').write_text(json.dumps(m,indent=2));print(json.dumps(m,indent=2))
