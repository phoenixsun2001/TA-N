from pathlib import Path
import json,warnings
import geopandas as g,shapely
from shapely.strtree import STRtree
from convert_gcj02 import convert
out=Path('outputs/notch_example_v9');out.mkdir(exist_ok=True)
d=g.read_file('outputs/closed_small_holes_v8/blocks_projected_final.gpkg').to_crs(32651);idx=d.index[d.piece_id==18274][0];original=d.geometry.iloc[idx];closed=original.buffer(8,quad_segs=8).buffer(-8,quad_segs=8).union(original);added=closed.difference(original);tree=STRtree(d.geometry.to_numpy());others=[i for i in tree.query(added,predicate='intersects') if i!=idx];overlap=sum(added.intersection(d.geometry.iloc[i]).area for i in others)
if overlap>.01:raise RuntimeError('Example overlaps adjacent blocks')
water=g.read_file('data/raw/osm/gis_osm_water_a_free_1.shp').to_crs(32651);wt=STRtree(water.geometry.to_numpy());wi=wt.query(added,predicate='intersects');water_area=sum(added.intersection(water.geometry.iloc[i]).area for i in wi)
if water_area>.01:raise RuntimeError('Example overlaps mapped water')
d.loc[idx,'geometry']=closed;d.loc[idx,'area_km2']=closed.area/1e6;d.loc[idx,'holes_remaining']=len(closed.interiors);d['notch_adjusted']=d.piece_id==18274
d.to_file(out/'blocks_projected.gpkg',driver='GPKG');g.GeoDataFrame({'piece_id':[18274],'geometry':[added]},crs=32651).to_file(out/'added_land_review.gpkg',driver='GPKG')
x=d.to_crs(4326);x.geometry=shapely.transform(x.geometry.array,convert,include_z=False);bad=~x.is_valid;x.loc[bad,'geometry']=shapely.make_valid(x.loc[bad,'geometry'].array);x=x.explode(index_parts=False,ignore_index=True);x=x[x.geom_type=='Polygon'];x.crs=None
with warnings.catch_warnings():
 warnings.filterwarnings('ignore',message="'crs' was not provided.*");x.to_file(out/'blocks_gcj02.geojson',driver='GeoJSON')
review=g.GeoDataFrame({'piece_id':[18274],'geometry':[added]},crs=32651).to_crs(4326);review.geometry=shapely.transform(review.geometry.array,convert);review.crs=None
with warnings.catch_warnings():
 warnings.filterwarnings('ignore',message="'crs' was not provided.*");review.to_file(out/'added_land_review_gcj02.geojson',driver='GeoJSON')
m={'piece_id':18274,'buffer_m':8,'area_added_m2':added.area,'adjacent_overlap_m2':overlap,'water_overlap_m2':water_area,'invalid':int((~d.is_valid).sum()),'gcj_invalid':int((~x.is_valid).sum()),'note':'Local trial closes the notch including a portion of mapped Xingfu Road. No citywide application.'};(out/'metrics.json').write_text(json.dumps(m,indent=2));print(json.dumps(m,indent=2))
