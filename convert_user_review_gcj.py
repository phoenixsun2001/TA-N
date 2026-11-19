from pathlib import Path
import json,warnings,shutil
import geopandas as g,shapely
from convert_gcj02 import convert
root=Path('outputs/user_data_test');out=root/'gcj02';out.mkdir(exist_ok=True)
items=[(Path('data/user_shanghai/national_roads_candidate_filtered.gpkg'),'user_road_candidates_gcj02.geojson'),(Path('data/user_shanghai/national_roads_water_review.gpkg'),'user_road_water_review_gcj02.geojson'),(Path('data/user_shanghai/rail_review_candidates.gpkg'),'user_rail_review_gcj02.geojson'),(Path('data/user_shanghai/hydrology_reference.gpkg'),'user_hydrology_reference_gcj02.geojson'),(Path('data/user_shanghai/pedestrian_review_candidates.gpkg'),'user_pedestrian_review_gcj02.geojson'),(root/'user_bridge_candidate_links.gpkg','user_bridge_candidate_links_gcj02.geojson'),(root/'repaired_passage_links.gpkg','passage_links_gcj02.geojson'),(root/'mall_anchor_audit.geojson','mall_anchor_audit_gcj02.geojson'),(Path('outputs/road_network_v2/inferred_gap_connectors.gpkg'),'inferred_gap_connectors_gcj02.geojson')];m=[]
for src,name in items:
 print('converting',src,flush=True);df=g.read_file(src).to_crs(4326);df.geometry=shapely.transform(df.geometry.array,convert,include_z=False);bad=~df.is_valid
 if bad.any():df.loc[bad,'geometry']=shapely.make_valid(df.loc[bad,'geometry'].array)
 df.crs=None
 with warnings.catch_warnings():
  warnings.filterwarnings('ignore',message="'crs' was not provided.*");df.to_file(out/name,driver='GeoJSON')
 m.append({'source':str(src),'output':name,'features':len(df),'validity_repairs':int(bad.sum())})
units=g.read_file(out/'road_graph_units_gcj02.geojson');pending=units[units.zone_type=='road_network_incomplete'].copy();pending.crs=None
with warnings.catch_warnings():
 warnings.filterwarnings('ignore',message="'crs' was not provided.*");pending.to_file(out/'unresolved_large_areas_gcj02.geojson',driver='GeoJSON')
shutil.copy2('outputs/road_network_v2/gcj02/commercial_core_candidates_gcj02.geojson',out/'commercial_core_candidates_gcj02.geojson');m.extend([{'source':'road_graph_units_gcj02.geojson','output':'unresolved_large_areas_gcj02.geojson','features':len(pending)},{'source':'previous_commercial_core','output':'commercial_core_candidates_gcj02.geojson','features':1981}]);(out/'support_manifest.json').write_text(json.dumps({'coordinates':'GCJ-02 approximate, not EPSG:4326','layers':m},ensure_ascii=False,indent=2));print(json.dumps(m,ensure_ascii=False,indent=2))
