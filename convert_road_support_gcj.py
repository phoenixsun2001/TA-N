from pathlib import Path
import json,warnings,shapely,geopandas as g
from convert_gcj02 import convert
root=Path('outputs/road_network_v2');out=root/'gcj02';out.mkdir(exist_ok=True)
items=[('commercial_core_candidates.gpkg','commercial_core_candidates_gcj02.geojson'),('inferred_gap_connectors.gpkg','inferred_gap_connectors_gcj02.geojson'),('repaired_passage_links.gpkg','passage_links_gcj02.geojson'),('mall_anchor_audit.geojson','mall_anchor_audit_gcj02.geojson')];summary=[]
for source,target in items:
 print('converting',source,flush=True);d=g.read_file(root/source).to_crs(4326);d.geometry=shapely.transform(d.geometry.array,convert,include_z=False);bad=~d.is_valid
 if bad.any():d.loc[bad,'geometry']=shapely.make_valid(d.loc[bad,'geometry'].array)
 d.crs=None
 with warnings.catch_warnings():
  warnings.filterwarnings('ignore',message="'crs' was not provided.*");d.to_file(out/target,driver='GeoJSON')
 summary.append({'source':source,'target':target,'features':len(d),'validity_repairs':int(bad.sum())})
# Inspect big unresolved areas independently; their boundary matches the shared-edge partition above.
units=g.read_file(out/'road_graph_units_gcj02.geojson');pending=units[units.zone_type=='road_network_incomplete'].copy();pending.crs=None
with warnings.catch_warnings():
 warnings.filterwarnings('ignore',message="'crs' was not provided.*");pending.to_file(out/'unresolved_large_areas_gcj02.geojson',driver='GeoJSON')
summary.append({'source':'shared_boundary_partition','target':'unresolved_large_areas_gcj02.geojson','features':len(pending)})
(out/'manifest.json').write_text(json.dumps({'crs':'GCJ-02 approximate','input_assumption':'WGS84; municipal boundary datum unverified','axis':'longitude,latitude','layers':summary},ensure_ascii=False,indent=2));print(json.dumps(summary,ensure_ascii=False,indent=2))
