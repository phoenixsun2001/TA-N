from pathlib import Path
import json,warnings
import numpy as np,pandas as pd,geopandas as g,shapely
from shapely.geometry import Point,LineString
from shapely.ops import nearest_points,linemerge,split
from shapely.strtree import STRtree
from convert_gcj02 import convert
out=Path('outputs/poi_anchored_v16');out.mkdir(exist_ok=True)
def osm(n):return g.read_file('data/raw/osm/gis_osm_'+n+'_free_1.shp').to_crs(32651)
city=g.read_file('data/raw/shanghai_boundary.geojson').to_crs(32651).geometry.union_all();shapely.prepare(city)
items=[]
for source,classes in [('transport',['ferry_terminal']),('transport_a',['ferry_terminal']),('traffic',['marina','slipway'])]:
 d=osm(source);d=d[d.fclass.isin(classes)].copy()
 for r in d.itertuples():items.append({'poi_id':str(r.osm_id),'source':source,'name':r.name or '', 'kind':'ferry_or_marina','geometry':r.geometry if r.geometry.geom_type=='Point' else r.geometry.representative_point()})
d=osm('traffic');d=d[(d.fclass=='parking')&d.name.fillna('').str.contains('入口|出口|出入口')]
for r in d.itertuples():items.append({'poi_id':str(r.osm_id),'source':'traffic','name':r.name or '','kind':'parking_entrance','geometry':r.geometry})
p=g.GeoDataFrame(items,crs=32651);p=p[shapely.intersects(city,p.geometry.array)].copy().reset_index(drop=True)
blocks=g.read_file('outputs/riverfront_v15/blocks_projected.gpkg').to_crs(32651);bg=blocks.geometry.to_numpy();bt=STRtree(bg)
roads=osm('roads');extra=g.read_file('data/user_shanghai/national_roads_candidate_filtered.gpkg').to_crs(32651);roads=g.GeoDataFrame(pd.concat([roads,extra],ignore_index=True),crs=32651);roads=roads[roads.fclass.isin(['motorway','motorway_link','trunk','trunk_link','primary','primary_link','secondary','secondary_link','tertiary','tertiary_link','residential','unclassified'])&(roads.tunnel!='T')&(roads.bridge!='T')&(roads.layer.fillna(0)==0)].reset_index(drop=True);rg=roads.geometry.to_numpy();rt=STRtree(rg)
by_block={};audits=[]
for ix,r in p.iterrows():
 near=bt.query(r.geometry,predicate='dwithin',distance=20)
 if not len(near):audits.append('no_block_within_20m');continue
 target=min(near,key=lambda i:bg[i].distance(r.geometry));by_block.setdefault(int(target),[]).append(ix);audits.append('queued')
p['status']=audits
rows=[];cutrows=[]
for i,row in blocks.iterrows():
 attrs=row.drop(labels='geometry').to_dict();geom=row.geometry;parts=[geom]
 for ix in by_block.get(i,[]):
  poi=p.geometry.iloc[ix];near=rt.query(geom,predicate='intersects');clipped=[z for j in near for z in shapely.get_parts(rg[j].intersection(geom)) if z.geom_type=='LineString' and z.length>1]
  try:merged=linemerge(shapely.union_all(clipped)) if len(clipped)>1 else clipped[0] if clipped else None
  except Exception:merged=None
  options=[]
  if merged is not None:
   for way in shapely.get_parts(merged):
    if way.geom_type!='LineString' or way.length<30 or way.distance(poi)>30:continue
    a,b=Point(way.coords[0]),Point(way.coords[-1]);da,db=geom.boundary.distance(a),geom.boundary.distance(b)
    if da>20 or db>20:continue
    aa=nearest_points(a,geom.boundary)[1];bb=nearest_points(b,geom.boundary)[1];coords=[(aa.x,aa.y),*way.coords,(bb.x,bb.y)];va=np.asarray(coords[0])-np.asarray(coords[1]);vb=np.asarray(coords[-1])-np.asarray(coords[-2]);na=np.linalg.norm(va);nb=np.linalg.norm(vb)
    if na>1e-6:coords[0]=tuple(np.asarray(coords[0])+va/na*.05)
    if nb>1e-6:coords[-1]=tuple(np.asarray(coords[-1])+vb/nb*.05)
    options.append((way.distance(poi),LineString(coords),max(da,db)))
  options.sort(key=lambda z:z[0]);found=False
  for dist,line,gap in options:
   new=[];changed=False
   for part in parts:
    try:pieces=[z for z in split(part,line).geoms if z.geom_type=='Polygon' and z.area>1]
    except Exception:pieces=[]
    if len(pieces)>1 and min(z.area for z in pieces)>=5000 and abs(sum(z.area for z in pieces)-part.area)<.1:new.extend(pieces);changed=True
    else:new.append(part)
   if changed:
    parts=new;p.loc[ix,'status']='split_candidate';cutrows.append({'poi_id':str(p.poi_id.iloc[ix]),'poi_kind':str(p.kind.iloc[ix]),'poi_name':str(p.name.iloc[ix]),'distance_m':dist,'max_connector_m':gap,'geometry':line.intersection(geom)});found=True;break
  if not found:p.loc[ix,'status']='no_qualified_road_cut'
 for part in parts:rows.append({**attrs,'poi_cut':len(parts)>1,'geometry':part})
r=g.GeoDataFrame(rows,crs=32651);r['poi_piece_id']=np.arange(1,len(r)+1);r['area_km2']=r.area/1e6;r.to_file(out/'blocks_projected.gpkg',driver='GPKG');p.to_file(out/'poi_audit.gpkg',driver='GPKG')
if cutrows:g.GeoDataFrame(cutrows,crs=32651).to_file(out/'poi_cut_lines.gpkg',driver='GPKG')
for frame,name in [(r,'blocks_gcj02.geojson'),(p,'poi_audit_gcj02.geojson')]+([(g.GeoDataFrame(cutrows,crs=32651),'poi_cut_lines_gcj02.geojson')] if cutrows else []):
 x=frame.to_crs(4326);x.geometry=shapely.transform(x.geometry.array,convert,include_z=False);bad=~x.is_valid;x.loc[bad,'geometry']=shapely.make_valid(x.loc[bad,'geometry'].array);x=x.explode(index_parts=False,ignore_index=True)
 if name.startswith('blocks'):x=x[x.geom_type=='Polygon'];x.geometry=[shapely.Polygon(z.exterior) for z in x.geometry]
 x.crs=None
 with warnings.catch_warnings():
  warnings.filterwarnings('ignore',message="'crs' was not provided.*");x.to_file(out/name,driver='GeoJSON')
 if name.startswith('blocks'):gcj={'features':len(x),'invalid':int((~x.is_valid).sum()),'holes':sum(len(z.interiors) for z in x.geometry)}
m={'poi_count':len(p),'poi_kind':p.kind.value_counts().to_dict(),'poi_status':p.status.value_counts().to_dict(),'cut_lines':len(cutrows),'input_blocks':len(blocks),'output_blocks':len(r),'area_delta_m2':float(r.area.sum()-blocks.area.sum()),'invalid':int((~r.is_valid).sum()),'holes':sum(len(z.interiors) for z in r.geometry),'gcj':gcj};(out/'metrics.json').write_text(json.dumps(m,ensure_ascii=False,indent=2));print(json.dumps(m,ensure_ascii=False,indent=2))
