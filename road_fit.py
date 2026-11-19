from pathlib import Path
import json, os
import numpy as np, geopandas as gpd, shapely
from shapely.geometry import Polygon
from shapely import STRtree

root=Path('outputs/city_imagery');out=root/os.getenv('ROAD_FIT_OUTPUT','road_fit');out.mkdir(exist_ok=True)
src=gpd.read_file(root/'closed/trade_areas_no_holes.geojson').to_crs(32651)
roads=gpd.read_file('data/raw/osm/gis_osm_roads_free_1.shp').to_crs(32651)
roads=roads[roads.fclass.isin(['motorway','trunk','primary','secondary'])&(roads.bridge!='T')&(roads.tunnel!='T')&(roads.layer.fillna(0)==0)].copy()
# Clip candidate road search to the neighbourhood of commercial areas.
roads=roads[roads.intersects(src.geometry.union_all().buffer(500))].copy(); roads=roads.explode(index_parts=False).reset_index(drop=True)
roads=roads[roads.geom_type=='LineString'].reset_index(drop=True);rg=roads.geometry.to_numpy();tree=STRtree(rg)
boundary=gpd.read_file('data/raw/shanghai_boundary.geojson').to_crs(32651).geometry.union_all().boundary
maxdist=float(os.getenv('ROAD_FIT_DISTANCE','150'));interval=50; cos_min=np.cos(np.deg2rad(float(os.getenv('ROAD_FIT_ANGLE','40')))); result=[]; changed=[]
for row in src.itertuples():
 geom=shapely.segmentize(row.geometry,interval); coords=np.asarray(geom.exterior.coords);points=shapely.points(coords[:-1]);prev=np.roll(coords[:-1],1,axis=0);nxt=np.roll(coords[:-1],-1,axis=0);tangent=nxt-prev;norm=np.linalg.norm(tangent,axis=1);tangent=np.divide(tangent,norm[:,None],out=np.zeros_like(tangent),where=norm[:,None]>0)
 pair,dist=tree.query_nearest(points,max_distance=maxdist,return_distance=True,all_matches=False);pi,ri=pair; lines=rg[ri]; loc=shapely.line_locate_point(lines,points[pi]);line_length=shapely.length(lines);a=shapely.line_interpolate_point(lines,np.maximum(loc-10,0));b=shapely.line_interpolate_point(lines,np.minimum(loc+10,line_length));vec=np.column_stack([shapely.get_x(b)-shapely.get_x(a),shapely.get_y(b)-shapely.get_y(a)]);vn=np.linalg.norm(vec,axis=1);vec=np.divide(vec,vn[:,None],out=np.zeros_like(vec),where=vn[:,None]>0);dot=np.abs(np.sum(vec*tangent[pi],axis=1));allowed=(dot>=cos_min)&(shapely.distance(points[pi],boundary)>1)
 chosen=pi[allowed];matched=shapely.line_interpolate_point(lines[allowed],loc[allowed]);coords[chosen,0]=shapely.get_x(matched);coords[chosen,1]=shapely.get_y(matched);coords[-1]=coords[0];p=Polygon(coords);result.append(p);changed.append(len(chosen))
r=src.copy();r.geometry=result;r['snapped_points']=changed
bad=~r.is_valid; valid=r[~bad]; union=valid.geometry.union_all();overlap=valid.area.sum()-union.area
print(json.dumps({'roads':len(roads),'points_snapped':int(sum(changed)),'features_changed':int(sum(np.array(changed)>0)),'invalid':int(bad.sum()),'overlap_m2':float(overlap)},indent=2),flush=True)
r.to_file(out/'trial_candidates.gpkg',driver='GPKG')
