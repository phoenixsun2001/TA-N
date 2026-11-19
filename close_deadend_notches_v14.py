from pathlib import Path
import json,warnings
import geopandas as g,numpy as np,shapely
from shapely.strtree import STRtree
from convert_gcj02 import convert
out=Path('outputs/deadend_notches_v14');out.mkdir(exist_ok=True)
d=g.read_file('outputs/fill_exterior_v13/blocks_projected.gpkg').to_crs(32651);geoms=d.geometry.to_numpy().copy();tree=STRtree(geoms)
water=g.read_file('data/raw/osm/gis_osm_water_a_free_1.shp').to_crs(32651);wg=shapely.make_valid(water.geometry.array);wt=STRtree(wg)
# Estimate the size of internal road notches from the input corridor width.
added=[];skips={};cells={}
def cell_keys(z):
 a,b,c,e=z.bounds
 return [(x,y) for x in range(int(np.floor(a/500)),int(np.floor(c/500))+1) for y in range(int(np.floor(b/500)),int(np.floor(e/500))+1)]
for i,p in enumerate(geoms):
 if p.area<2000:continue
 try:diff=p.buffer(10,quad_segs=4).buffer(-10,quad_segs=4).union(p).difference(p)
 except Exception:skips['geometry_error']=skips.get('geometry_error',0)+1;continue
 for z in shapely.get_parts(diff):
  if z.geom_type!='Polygon' or z.area<50:continue
  minx,miny,maxx,maxy=z.bounds;long=max(maxx-minx,maxy-miny);short=min(maxx-minx,maxy-miny)
  if long<20 or long<short*1.5:continue
  near=tree.query(z,predicate='intersects');other=sum(z.intersection(geoms[j]).area for j in near if j!=i)
  if other>.01:skips['other_block']=skips.get('other_block',0)+1;continue
  water_overlap=sum(z.intersection(wg[j]).area for j in wt.query(z,predicate='intersects'))
  if water_overlap>.01:skips['water']=skips.get('water',0)+1;continue
  keys=cell_keys(z);near_added=set(j for key in keys for j in cells.get(key,[]))
  if any(added[j]['geometry'].intersects(z) and added[j]['geometry'].intersection(z).area>.01 for j in near_added if added[j]['block_id']!=int(d.block_id.iloc[i])):skips['prior_added']=skips.get('prior_added',0)+1;continue
  joined=geoms[i].union(z)
  if joined.geom_type!='Polygon' or not joined.is_valid or len(joined.interiors):skips['geometry']=skips.get('geometry',0)+1;continue
  delta=joined.difference(geoms[i]);geoms[i]=joined
  if delta.area>.001:
   index=len(added);added.append({'block_id':int(d.block_id.iloc[i]),'area_m2':delta.area,'geometry':delta})
   for key in keys:cells.setdefault(key,[]).append(index)
d.geometry=geoms;d['area_km2']=d.area/1e6;d['notch_closed']=d.block_id.isin([r['block_id'] for r in added]);d.to_file(out/'blocks_projected.gpkg',driver='GPKG')
if added:g.GeoDataFrame(added,crs=32651).to_file(out/'added_areas.gpkg',driver='GPKG')
for frame,name in [(d,'blocks_gcj02.geojson'),(g.GeoDataFrame(added,crs=32651),'added_areas_gcj02.geojson')]:
 if not len(frame):continue
 x=frame.to_crs(4326);x.geometry=shapely.transform(x.geometry.array,convert,include_z=False);bad=~x.is_valid;x.loc[bad,'geometry']=shapely.make_valid(x.loc[bad,'geometry'].array);x=x.explode(index_parts=False,ignore_index=True);x=x[x.geom_type=='Polygon'];x.crs=None
 with warnings.catch_warnings():
  warnings.filterwarnings('ignore',message="'crs' was not provided.*");x.to_file(out/name,driver='GeoJSON')
 if name.startswith('blocks'):gx={'features':len(x),'invalid':int((~x.is_valid).sum()),'holes':sum(len(p.interiors) for p in x.geometry)}
m={'source_features':len(d),'closed_candidates':len(added),'modified_blocks':int(d.notch_closed.sum()),'added_area_km2':sum(x['area_m2'] for x in added)/1e6,'skipped':skips,'invalid':int((~d.is_valid).sum()),'holes':sum(len(p.interiors) for p in d.geometry),'gcj':gx};(out/'metrics.json').write_text(json.dumps(m,indent=2));print(json.dumps(m,indent=2))
