from pathlib import Path
import json,warnings
import numpy as np, geopandas as g,shapely
from shapely.geometry import Point,LineString,Polygon
from shapely.ops import nearest_points,split
from shapely.strtree import STRtree
from convert_gcj02 import convert
out=Path('outputs/poi_direct_v17');out.mkdir(exist_ok=True)
d=g.read_file('outputs/riverfront_v15/blocks_projected.gpkg').to_crs(32651);pois=g.read_file('outputs/poi_anchored_v16/poi_audit_extended.gpkg').to_crs(32651);geo=d.geometry.to_numpy();tree=STRtree(geo);groups={};statuses=['']*len(pois)
for i,p in enumerate(pois.geometry):
 matches=tree.query(p,predicate='dwithin',distance=200)
 if not len(matches):statuses[i]='no_land_within_200m';continue
 suitable=[int(j) for j in matches if geo[j].area>=1000]
 index=min(suitable or [int(j) for j in matches],key=lambda j:geo[j].distance(p));groups.setdefault(index,[]).append(i)
rows=[];cuts=[]
for i,row in d.iterrows():
 attrs=row.drop(labels='geometry').to_dict();parts=[row.geometry]
 for pi in groups.get(i,[]):
  point=pois.geometry.iloc[pi];target=min(range(len(parts)),key=lambda n:parts[n].distance(point));p=parts[target];nearest=nearest_points(p,point)[0];rep=p.representative_point();dist=p.distance(point);made=False
  # Directions favour a cross-section of the long axis, then explore alternatives.
  box=p.minimum_rotated_rectangle;xy=np.asarray(box.exterior.coords);vec=np.diff(xy,axis=0);lens=np.linalg.norm(vec,axis=1);base=np.arctan2(vec[np.argmin(lens),1],vec[np.argmin(lens),0]);span=max(p.bounds[2]-p.bounds[0],p.bounds[3]-p.bounds[1])*3+1
  for offset in [0,.5,1,2,5,10,20,50,100,200]:
   center=point if p.contains(point) else Point(nearest.x+(rep.x-nearest.x)*min(1,offset/max(nearest.distance(rep),1e-9)),nearest.y+(rep.y-nearest.y)*min(1,offset/max(nearest.distance(rep),1e-9)))
   if not p.covers(center):continue
   for angle in base+np.linspace(0,np.pi,24,endpoint=False):
    ux=np.cos(angle);uy=np.sin(angle);line=LineString([(center.x-ux*span,center.y-uy*span),(center.x+ux*span,center.y+uy*span)])
    try:q=[z for z in split(p,line).geoms if z.geom_type=='Polygon' and z.area>1]
    except Exception:continue
    if len(q)>1 and abs(sum(z.area for z in q)-p.area)<.1 and all(z.is_valid and not z.interiors for z in q):
     parts[target:target+1]=q;cuts.append({'poi_id':str(pois.poi_id.iloc[pi]),'poi_kind':str(pois.kind.iloc[pi]),'name':str(pois.name.iloc[pi]),'poi_to_block_m':float(dist),'poi_to_cut_m':float(point.distance(line)),'geometry':line.intersection(p)});statuses[pi]='split';made=True;break
   if made:break
  if not made:statuses[pi]='no_valid_cut'
 for p in parts:rows.append({**attrs,'poi_direct_cut':len(parts)>1,'geometry':p})
r=g.GeoDataFrame(rows,crs=32651);r['direct_piece_id']=np.arange(1,len(r)+1);r['area_km2']=r.area/1e6;r.to_file(out/'blocks_projected.gpkg',driver='GPKG');pois['direct_status']=statuses;pois.to_file(out/'poi_audit.gpkg',driver='GPKG')
cut=g.GeoDataFrame(cuts,crs=32651);cut.to_file(out/'cut_lines.gpkg',driver='GPKG')
for frame,name in [(r,'blocks_gcj02.geojson'),(pois,'poi_audit_gcj02.geojson'),(cut,'cut_lines_gcj02.geojson')]:
 x=frame.to_crs(4326);x.geometry=shapely.transform(x.geometry.array,convert,include_z=False);bad=~x.is_valid;x.loc[bad,'geometry']=shapely.make_valid(x.loc[bad,'geometry'].array);x=x.explode(index_parts=False,ignore_index=True)
 if name.startswith('blocks'):x=x[x.geom_type=='Polygon'];x.geometry=[Polygon(z.exterior) for z in x.geometry]
 x.crs=None
 with warnings.catch_warnings():
  warnings.filterwarnings('ignore',message="'crs' was not provided.*");x.to_file(out/name,driver='GeoJSON')
 if name.startswith('blocks'):gcj={'features':len(x),'invalid':int((~x.is_valid).sum()),'holes':sum(len(z.interiors) for z in x.geometry)}
m={'poi':len(pois),'status':pois.direct_status.value_counts().to_dict(),'cut_lines':len(cuts),'input_blocks':len(d),'output_blocks':len(r),'area_delta_m2':float(r.area.sum()-d.area.sum()),'invalid':int((~r.is_valid).sum()),'holes':sum(len(z.interiors) for z in r.geometry),'max_poi_cut_distance_m':float(max([z['poi_to_cut_m'] for z in cuts],default=0)),'gcj':gcj};(out/'metrics.json').write_text(json.dumps(m,ensure_ascii=False,indent=2));print(json.dumps(m,ensure_ascii=False,indent=2))
