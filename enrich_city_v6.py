from pathlib import Path
import json,warnings
import geopandas as g
import numpy as np
import shapely,rasterio
from shapely.strtree import STRtree
from shapely.geometry import Point,LineString
from convert_gcj02 import convert
out=Path('outputs/rebuild_v6')
d=g.read_file(out/'land_blocks_projected.gpkg').to_crs(32651)
d=d[['block_id','geometry']].copy();d['area_km2']=d.area/1e6
points=g.read_file('data/raw/osm/gis_osm_pois_free_1.shp').to_crs(32651)
join=g.sjoin(points[['geometry']],d[['block_id','geometry']],predicate='within');counts=join.groupby('block_id').size();d['poi_count']=d.block_id.map(counts).fillna(0).astype(int)
mall=g.read_file('data/raw/osm/gis_osm_pois_a_free_1.shp').to_crs(32651);mall=mall[mall.fclass=='mall'].copy();mall.geometry=mall.geometry.representative_point();joined=g.sjoin(mall[['geometry']],d[['block_id','geometry']],predicate='within');counts=joined.groupby('block_id').size();d['mall_count']=d.block_id.map(counts).fillna(0).astype(int)
with rasterio.open('outputs/city_imagery/city_features.tif') as ds:
 reps=d.geometry.representative_point();bands=['viirs_index','ndvi','ndbi','s1_vv_db'];ids=[ds.descriptions.index(n)+1 for n in bands];samples=np.array(list(ds.sample(zip(reps.x,reps.y),indexes=ids)))
 for i,n in enumerate(bands):d[n]=np.where((samples[:,i]==-9999)|~np.isfinite(samples[:,i]),np.nan,samples[:,i])
d['review_small']=d.area<1000;d['review_large']=d.area>2e6;d['has_holes']=[bool(len(p.interiors)) for p in d.geometry];d['classification']='unclassified';d.loc[d.poi_count>=5,'classification']='poi_concentration_candidate';d.loc[d.mall_count>0,'classification']='mall_candidate';d['imagery_sample']='representative_point_250m'
d.to_file(out/'reconstructed_blocks.gpkg',driver='GPKG')
x=d.to_crs(4326);x.geometry=shapely.transform(x.geometry.array,convert,include_z=False);bad=~x.is_valid;x.loc[bad,'geometry']=shapely.make_valid(x.loc[bad,'geometry'].array);x=x.explode(index_parts=False,ignore_index=True);x=x[x.geom_type=='Polygon'].copy();x.crs=None
with warnings.catch_warnings():
 warnings.filterwarnings('ignore',message="'crs' was not provided.*");x.to_file(out/'reconstructed_blocks_gcj02.geojson',driver='GeoJSON')
# Existing explicitly mapped bridge links are candidates. Endpoints are snapped to blocks; they are not a verified pedestrian routing network.
links=g.read_file('outputs/user_data_test/repaired_passage_links.gpkg').to_crs(32651);links=links[links.kind=='bridge_link'];geoms=d.geometry.to_numpy();tree=STRtree(geoms);rows=[]
for r in links.itertuples():
 z=r.geometry
 if z.geom_type!='LineString':continue
 ends=[Point(z.coords[0]),Point(z.coords[-1])];matches=[]
 for p in ends:
  near=tree.query(p,predicate='dwithin',distance=35)
  if not len(near):break
  matches.append(int(near[np.argmin(shapely.distance(geoms[near],p))]))
 if len(matches)==2 and matches[0]!=matches[1]:rows.append({'block_a':int(d.block_id.iloc[matches[0]]),'block_b':int(d.block_id.iloc[matches[1]]),'status':'mapped_crossing_unverified','geometry':z})
if rows:
 l=g.GeoDataFrame(rows,crs=32651);l.to_file(out/'crossing_candidates.gpkg',driver='GPKG');l=l.to_crs(4326);l.geometry=shapely.transform(l.geometry.array,convert,include_z=False);l.crs=None
 with warnings.catch_warnings():
  warnings.filterwarnings('ignore',message="'crs' was not provided.*");l.to_file(out/'crossing_candidates_gcj02.geojson',driver='GeoJSON')
m={'blocks':len(d),'valid':int(d.is_valid.sum()),'holes':int(d.has_holes.sum()),'small_under_1000m2':int(d.review_small.sum()),'large_over_2km2':int(d.review_large.sum()),'large_area_km2':float(d.area[d.review_large].sum()/1e6),'mall_points_assigned':len(joined),'mapped_crossing_candidates':len(rows),'gcj_features':len(x),'gcj_invalid':int((~x.is_valid).sum())};(out/'enrichment_metrics.json').write_text(json.dumps(m,indent=2));print(json.dumps(m,indent=2))
