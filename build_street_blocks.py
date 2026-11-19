from pathlib import Path
import json,time
import geopandas as gpd, numpy as np, shapely
from shapely.geometry import box,Polygon
from shapely.ops import polygonize,unary_union
from pyproj import Transformer
out=Path('outputs/block_fit');out.mkdir(parents=True,exist_ok=True)
def read(name):return gpd.read_file(f'data/raw/osm/gis_osm_{name}_free_1.shp').to_crs(32651)
roads=read('roads');rail=read('railways');water=read('water_a');traffic=read('traffic');centers=gpd.read_file('outputs/centers.geojson').to_crs(32651)
# Cut lines are physical ground-level major and neighbourhood roads. Grade-separated bridges/tunnels are excluded.
cut=roads[roads.fclass.isin(['motorway','trunk','primary','secondary','tertiary'])&(roads.bridge!='T')&(roads.tunnel!='T')&(roads.layer.fillna(0)==0)]
major=cut[cut.fclass.isin(['motorway','trunk','primary','secondary'])]
rail=rail[rail.fclass.isin(['rail','light_rail','tram'])&(rail.bridge!='T')&(rail.tunnel!='T')&(rail.layer.fillna(0)==0)]
water=water[water.fclass.isin(['water','riverbank','reservoir','dock'])]
walk=roads[roads.fclass.isin(['footway','path','pedestrian','steps'])]
bridges=walk[(walk.bridge=='T')|(walk.tunnel=='T')]
cross=traffic[traffic.fclass.isin(['pedestrian_crossing','railway_crossing'])]
projection=Transformer.from_crs(4326,32651,always_xy=True)
sites={'xujiahui':(121.4365,31.195),'lujiazui':(121.505,31.237),'jinqiao':(121.603,31.252)}
for name,coord in sites.items():
 x,y=projection.transform(*coord);domain=box(x-2000,y-2000,x+2000,y+2000);bounds=domain.bounds
 def subset(d):return d.cx[bounds[0]:bounds[2],bounds[1]:bounds[3]]
 localcut=subset(cut);localrail=subset(rail);localwater=subset(water);localmajor=subset(major);localcross=subset(cross);localwalk=subset(walk);localbridges=subset(bridges)
 watergeom=unary_union([g.intersection(domain) for g in localwater.geometry if g.intersects(domain)]);watergeom=shapely.make_valid(watergeom)
 lines=[domain.boundary]+[g.intersection(domain) for g in localcut.geometry if g.intersects(domain)]+[g.intersection(domain) for g in localrail.geometry if g.intersects(domain)]
 if not watergeom.is_empty:lines.append(watergeom.boundary)
 print(name,'noding',len(lines),flush=True);network=unary_union(lines);faces=[p for p in polygonize(network) if p.area>0.05 and domain.covers(p.representative_point())]
 blocks=gpd.GeoDataFrame({'id':range(len(faces))},geometry=faces,crs=32651);blocks['is_water']=[not watergeom.is_empty and watergeom.covers(g.representative_point()) for g in faces]
 print(name,'blocks',len(blocks),'water',int(blocks.is_water.sum()),'valid',int((~blocks.is_valid).sum()),'coverage_gap',domain.difference(blocks.geometry.union_all()).area,flush=True)
 blocks.to_file(out/f'{name}_blocks.gpkg',driver='GPKG')
 data={'name':name,'domain_bounds':bounds,'cuts':len(localcut),'major_roads':len(localmajor),'rail':len(localrail),'water':len(localwater),'crossings':len(localcross),'pedestrian_ways':len(localwalk),'pedestrian_bridges_or_tunnels':len(localbridges),'blocks':len(blocks),'water_blocks':int(blocks.is_water.sum()),'centers':len(subset(centers))}
 (out/f'{name}_input_metrics.json').write_text(json.dumps(data,indent=2));print(json.dumps(data),flush=True)
