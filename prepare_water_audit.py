from pathlib import Path
import json,geopandas as g,numpy as np,shapely,pandas as pd
from shapely.strtree import STRtree
out=Path('data/user_shanghai');city=g.read_file('data/raw/shanghai_boundary.geojson').to_crs(32651).geometry.union_all();scope=city.buffer(2000);lines=g.read_file('data/raw/osm/gis_osm_waterways_free_1.shp').to_crs(32651);polygons=g.read_file('data/raw/osm/gis_osm_water_a_free_1.shp').to_crs(32651);tree=STRtree(np.concatenate([lines.geometry.to_numpy(),polygons.geometry.to_numpy()]));rows=[];metrics=[]
for fp in sorted(out.glob('我国*河流.shp.gpkg')):
 d=g.read_file(fp).to_crs(32651);geoms=shapely.intersection(d.geometry.array,scope);length=shapely.length(geoms);samples=[]
 for i in np.flatnonzero(length>50):
  n=max(1,int(np.ceil(length[i]/250)));samples.extend(shapely.line_interpolate_point(geoms[i],(np.arange(n)+.5)*length[i]/n))
 if samples:
  _,dist=tree.query_nearest(np.array(samples,dtype=object),return_distance=True,all_matches=False);metrics.append({'source':fp.name,'sample_points':len(samples),'spacing_m':250,'nearest_water_median_m':float(np.median(dist)),'nearest_water_p90_m':float(np.quantile(dist,.9))})
 for i in np.flatnonzero(length>0):rows.append({'source':fp.name,'name':d.NAME.iloc[i],'geometry':geoms[i]})
f=g.GeoDataFrame(rows,crs=32651);f.to_file(out/'hydrology_reference.gpkg',driver='GPKG');(out/'hydrology_audit.json').write_text(json.dumps(metrics,ensure_ascii=False,indent=2));print(json.dumps(metrics,ensure_ascii=False,indent=2))
