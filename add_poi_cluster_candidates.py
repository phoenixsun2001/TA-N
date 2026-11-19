from pathlib import Path
import geopandas as g,numpy as np,shapely,pandas as pd,json
from sklearn.cluster import DBSCAN
from shapely.ops import unary_union
out=Path('outputs/revised_poc');p=g.read_file('outputs/commercial_poi.geojson').to_crs(32651);base=g.read_file(out/'commercial_core_candidates.geojson').to_crs(32651);base=base[base.source_type!='osm_poi_cluster'].copy();u=base.geometry.union_all();shapely.prepare(u);outside=p[~shapely.covers(u,p.geometry.array)].copy();xy=np.column_stack([outside.geometry.x,outside.geometry.y]);outside['cluster']=DBSCAN(eps=150,min_samples=4).fit_predict(xy)
water=g.read_file('data/raw/osm/gis_osm_water_a_free_1.shp').to_crs(32651);water=water[water.fclass.isin(['water','riverbank','reservoir','dock'])].geometry.union_all();new=[]
for id,grp in outside[outside.cluster>=0].groupby('cluster'):
 geom=unary_union(grp.geometry.buffer(55)).difference(water);geom=shapely.make_valid(geom)
 if geom.is_empty:continue
 new.append({'source_id':f'poi_cluster_{id}','name':'','source_type':'osm_poi_cluster','poi_count':len(grp),'geometry':geom})
extra=g.GeoDataFrame(new,crs=32651);merged=g.GeoDataFrame(pd.concat([base,extra],ignore_index=True),crs=32651);merged.to_crs(4326).to_file(out/'commercial_core_candidates.geojson',driver='GeoJSON');extra.to_crs(4326).to_file(out/'additional_poi_clusters.geojson',driver='GeoJSON');metrics={'outside_poi':len(outside),'additional_clusters':len(extra),'poi_in_additional_clusters':int(extra.poi_count.sum()),'outside_poi_noise':int((outside.cluster==-1).sum()),'dbscan_eps_m':150,'min_samples':4,'buffer_m':55};(out/'poi_cluster_metrics.json').write_text(json.dumps(metrics,indent=2));print(json.dumps(metrics,indent=2))
