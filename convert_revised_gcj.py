from pathlib import Path
import geopandas as g,shapely,json,warnings
from convert_gcj02 import convert
root=Path('outputs/revised_poc');out=root/'gcj02';out.mkdir(exist_ok=True)
items=[('commercial_core_candidates.geojson','commercial_core_candidates_gcj02.geojson'),('additional_poi_clusters.geojson','additional_poi_clusters_gcj02.geojson'),('activity_units_parent.gpkg','activity_units_parent_gcj02.geojson'),('activity_partition.gpkg','activity_partition_gcj02.geojson'),('activity_components_no_holes.gpkg','activity_components_no_holes_gcj02.geojson'),('inferred_crossing_candidates.gpkg','inferred_crossings_gcj02.geojson'),('mall_anchor_audit.geojson','mall_anchor_audit_gcj02.geojson')];manifest=[]
for source,target in items:
 print('converting',source,flush=True);df=g.read_file(root/source).to_crs(4326);df.geometry=shapely.transform(df.geometry.array,convert,include_z=False);bad=~df.is_valid
 # Coordinate transformation can collapse microscopic artifacts; repair only affected geometries.
 if bad.any():df.loc[bad,'geometry']=shapely.make_valid(df.loc[bad,'geometry'].array)
 assert df.is_valid.all();df.crs=None
 with warnings.catch_warnings():
  warnings.filterwarnings('ignore',message="'crs' was not provided.*");df.to_file(out/target,driver='GeoJSON')
 manifest.append({'source':source,'target':target,'features':len(df),'validity_repairs':int(bad.sum())})
(out/'manifest.json').write_text(json.dumps({'coordinate_system':'GCJ-02 approximate','axis_order':'longitude latitude','source_assumption':'WGS84 assumed; boundary datum unverified','layers':manifest},ensure_ascii=False,indent=2));print(json.dumps(manifest,ensure_ascii=False,indent=2),flush=True)
