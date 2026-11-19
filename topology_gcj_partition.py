from pathlib import Path
import json, warnings
import geopandas as g,numpy as np,shapely
from shapely.ops import polygonize
from shapely.strtree import STRtree
from pyproj import Transformer
from convert_gcj02 import convert
root=Path('outputs/revised_poc');out=root/'gcj02';src=g.read_file(root/'activity_partition.gpkg').to_crs(32651);source_geom=src.geometry.to_numpy();print('noding shared boundaries',len(src),flush=True);network=shapely.union_all(shapely.boundary(source_geom));print('projecting network',flush=True);wg=g.GeoSeries([network],crs=32651).to_crs(4326).iloc[0];gn=shapely.transform(wg,convert);print('polygonizing',flush=True);faces=list(polygonize(gn));print('faces',len(faces),flush=True)
points=np.array([[v.representative_point().x,v.representative_point().y] for v in faces]);wgs=points.copy()
for _ in range(8):wgs-=convert(wgs)-points
trans=Transformer.from_crs(4326,32651,always_xy=True);x,y=trans.transform(wgs[:,0],wgs[:,1]);orig_points=shapely.points(x,y);tree=STRtree(source_geom);matches=tree.query(orig_points,predicate='intersects');owners=np.full(len(faces),-1,dtype='int32')
for pi,si in matches.T:
 if owners[pi]<0:owners[pi]=si
unmatched=np.flatnonzero(owners<0)
if len(unmatched):
 near,dist=tree.query_nearest(orig_points[unmatched],max_distance=1,return_distance=True,all_matches=False)
 for pos,si in near.T:owners[unmatched[pos]]=si
unmatched=np.flatnonzero(owners<0);print('unmatched faces',len(unmatched),flush=True)
if len(unmatched):raise RuntimeError('Cannot assign transformed faces safely')
# Recombine polygonized faces into their original feature, keeping all shared edges from the same transformed network.
facearray=np.asarray(faces,dtype=object);order=np.argsort(owners,kind='stable');sorted_ids=owners[order];starts=np.r_[0,np.flatnonzero(np.diff(sorted_ids))+1];ends=np.r_[starts[1:],len(order)];attributes=src.drop(columns='geometry').to_dict('records');rows=[]
for start,end in zip(starts,ends):
 owner=int(sorted_ids[start]);geom=shapely.union_all(facearray[order[start:end]]);rows.append({**attributes[owner],'geometry':geom})
final=g.GeoDataFrame(rows,crs=None);print('features',len(final),'invalid',int((~final.is_valid).sum()),flush=True)
if not final.is_valid.all():raise RuntimeError('Invalid transformed partition')
with warnings.catch_warnings():
 warnings.filterwarnings('ignore',message="'crs' was not provided.*");final.to_file(out/'activity_partition_gcj02.geojson',driver='GeoJSON');final[final.unit_seed>0].to_file(out/'activity_units_parent_gcj02.geojson',driver='GeoJSON')
# A shared-network polygonization has no overlaps by construction; verify after reload in a metric pseudo-projection.
check=g.read_file(out/'activity_partition_gcj02.geojson').to_crs(32651);union=check.geometry.union_all();overlap=check.area.sum()-union.area;invalid=int((~check.is_valid).sum());print('overlap_m2',overlap,'invalid',invalid,flush=True)
if abs(overlap)>.1 or invalid:raise RuntimeError('Converted partition failed topology check')
metrics={'source_features':len(src),'faces':len(faces),'output_features':len(final),'original_features_absent':len(src)-len(final),'unmatched_faces':len(unmatched),'overlap_m2':float(overlap),'invalid_geometries':invalid,'note':'GCJ coordinates are treated as pseudo-WGS84 only for approximate metric QA.'};(out/'partition_topology_metrics.json').write_text(json.dumps(metrics,ensure_ascii=False,indent=2));print(json.dumps(metrics,indent=2),flush=True)
