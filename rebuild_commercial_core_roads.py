from pathlib import Path
import json
import geopandas as g,numpy as np,shapely,pandas as pd
from sklearn.cluster import DBSCAN
from shapely.ops import unary_union
out=Path('outputs/road_network_v2');p=g.read_file('outputs/commercial_poi.geojson').to_crs(32651);city=g.read_file('data/raw/shanghai_boundary.geojson').to_crs(32651).geometry.union_all();mall=g.read_file('data/raw/osm/gis_osm_pois_a_free_1.shp').to_crs(32651);mall=mall[(mall.fclass=='mall')&mall.geometry.representative_point().within(city)].copy();mall.geometry=shapely.make_valid(mall.geometry.array);mall['source_type']='osm_mall_aoi';mall['source_id']=mall.osm_id.astype(str);mall['poi_count']=-1;mall=mall[['source_id','name','source_type','poi_count','geometry']]
xy=np.column_stack((p.geometry.x,p.geometry.y));labels=DBSCAN(eps=150,min_samples=4).fit_predict(xy);water=g.read_file('data/raw/osm/gis_osm_water_a_free_1.shp').to_crs(32651);water=water[water.fclass.isin(['water','riverbank','reservoir','dock'])].geometry.union_all();poi_old=g.read_file('outputs/city_imagery/city_commercial_core.geojson').to_crs(32651).geometry.union_all();shapely.prepare(poi_old)
records=[]
for cluster in sorted(set(labels)-{-1}):
 sub=p.iloc[np.flatnonzero(labels==cluster)];ge=unary_union(sub.geometry.buffer(55)).intersection(city).difference(water);ge=shapely.make_valid(ge)
 if not ge.is_empty:records.append({'source_id':f'poi-cluster-{cluster}','name':'','source_type':'poi_cluster_150m','poi_count':len(sub),'geometry':ge})
noise=p.iloc[np.flatnonzero(labels==-1)].copy();noise=noise[shapely.covers(poi_old,noise.geometry.array)]
# Preserve points previously supported by the density surface without exporting its square pixel boundaries.
for row in noise.itertuples():
 geom=row.geometry.buffer(55).intersection(city).difference(water)
 if not geom.is_empty:records.append({'source_id':'density-poi-'+str(row.osm_id),'name':row.name or '','source_type':'density_supported_poi','poi_count':1,'geometry':geom})
cluster=g.GeoDataFrame(records,crs=32651);merged=g.GeoDataFrame(pd.concat([mall,cluster],ignore_index=True),crs=32651);merged.to_file(out/'commercial_core_candidates.gpkg',driver='GPKG');merged.to_crs(4326).to_file(out/'commercial_core_candidates.geojson',driver='GeoJSON');u=merged.geometry.union_all();shapely.prepare(u);metrics={'mall_aoi':len(mall),'poi_clusters':len(set(labels)-{-1}),'density_supported_noise':len(noise),'output_features':len(merged),'commercial_poi_total':len(p),'commercial_poi_covered':int(shapely.covers(u,p.geometry.array).sum()),'source_note':'POI buffer candidates are not verified trade area polygons; no raster grid is exported as their boundary.'};(out/'commercial_core_metrics.json').write_text(json.dumps(metrics,indent=2));print(json.dumps(metrics,indent=2))
