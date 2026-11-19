from pathlib import Path
import json,warnings
import geopandas as g,numpy as np,shapely
from shapely.geometry import Polygon
from shapely.strtree import STRtree
from convert_gcj02 import convert
out=Path('outputs/fill_exterior_v13');out.mkdir(exist_ok=True)
d=g.read_file('outputs/road_level_v12/reconstructed_blocks.gpkg').to_crs(32651);old=d.geometry.to_numpy();holes=sum(len(p.interiors) for p in old);d['original_holes']=[len(p.interiors) for p in old];new=np.array([Polygon(p.exterior) for p in old],dtype=object);d.geometry=new;d['added_area_m2']=shapely.area(new)-shapely.area(old);d['area_km2']=d.area/1e6;d.to_file(out/'blocks_projected.gpkg',driver='GPKG')
tree=STRtree(new);pairs=tree.query(new,predicate='intersects').T;pairs=pairs[pairs[:,0]<pairs[:,1]];intersections=shapely.area(shapely.intersection(new[pairs[:,0]],new[pairs[:,1]]));mask=intersections>.01
x=d.to_crs(4326);x.geometry=shapely.transform(x.geometry.array,convert,include_z=False);bad=~x.is_valid;x.loc[bad,'geometry']=shapely.make_valid(x.loc[bad,'geometry'].array);x=x.explode(index_parts=False,ignore_index=True);x=x[x.geom_type=='Polygon'];x.crs=None
with warnings.catch_warnings():
 warnings.filterwarnings('ignore',message="'crs' was not provided.*");x.to_file(out/'blocks_gcj02.geojson',driver='GeoJSON')
m={'source_features':len(d),'source_features_with_holes':int((d.original_holes>0).sum()),'source_hole_count':holes,'filled_area_km2':float(d.added_area_m2.sum()/1e6),'projected_invalid':int((~d.is_valid).sum()),'projected_holes':int(sum(len(p.interiors) for p in d.geometry)),'overlapping_pairs_over_0_01m2':int(mask.sum()),'overlap_pair_area_km2':float(intersections[mask].sum()/1e6),'gcj_features':len(x),'gcj_invalid':int((~x.is_valid).sum()),'gcj_holes':int(sum(len(p.interiors) for p in x.geometry))};(out/'metrics.json').write_text(json.dumps(m,indent=2));print(json.dumps(m,indent=2))
