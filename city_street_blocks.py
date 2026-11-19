from pathlib import Path
import json,time
import geopandas as gpd,shapely
from shapely.ops import polygonize,unary_union
out=Path('outputs/block_fit');out.mkdir(exist_ok=True)
def read(name):return gpd.read_file(f'data/raw/osm/gis_osm_{name}_free_1.shp').to_crs(32651)
city=gpd.read_file('data/raw/shanghai_boundary.geojson').to_crs(32651).geometry.union_all();bounds=city.bounds
roads=read('roads');rail=read('railways');water=read('water_a')
cut=roads[roads.fclass.isin(['motorway','trunk','primary','secondary','tertiary','residential','unclassified','service','living_street'])&(roads.bridge!='T')&(roads.tunnel!='T')&(roads.layer.fillna(0)==0)].cx[bounds[0]:bounds[2],bounds[1]:bounds[3]]
rail=rail[rail.fclass.isin(['rail','light_rail','tram'])&(rail.bridge!='T')&(rail.tunnel!='T')&(rail.layer.fillna(0)==0)].cx[bounds[0]:bounds[2],bounds[1]:bounds[3]]
water=water[water.fclass.isin(['water','riverbank','reservoir','dock'])].cx[bounds[0]:bounds[2],bounds[1]:bounds[3]]
print('source',len(cut),len(rail),len(water),flush=True)
watergeom=shapely.make_valid(water.geometry.union_all()); lines=list(cut.geometry)+list(rail.geometry)+[watergeom.boundary,city.boundary]
print('noding',len(lines),flush=True);network=unary_union(lines);print('polygonizing',flush=True);shapely.prepare(city)
faces=[p for p in polygonize(network) if p.area>0.05 and city.covers(p.representative_point())];print('faces',len(faces),flush=True)
shapely.prepare(watergeom);blocks=gpd.GeoDataFrame({'id':range(len(faces))},geometry=faces,crs=32651);blocks['is_water']=[watergeom.covers(p.representative_point()) for p in faces]
union=blocks.geometry.union_all();metrics={'road_lines':len(cut),'rail_lines':len(rail),'water_polygons':len(water),'faces':len(faces),'water_faces':int(blocks.is_water.sum()),'gap_m2':float(city.difference(union).area),'overlap_m2':float(blocks.area.sum()-union.area),'invalid':int((~blocks.is_valid).sum())};print(json.dumps(metrics,indent=2),flush=True)
blocks.to_file(out/'city_blocks.gpkg',driver='GPKG');(out/'city_input_metrics.json').write_text(json.dumps(metrics,indent=2))
