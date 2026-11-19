from pathlib import Path
import warnings
import numpy as np, geopandas as g, shapely
p=Path('outputs/barrier_cut_v3/land_blocks_gcj02.geojson');d=g.read_file(p);bad=~d.is_valid;print('invalid before',int(bad.sum()),flush=True)
d.loc[bad,'geometry']=shapely.make_valid(d.loc[bad,'geometry'].array)
d=d.explode(index_parts=False,ignore_index=True);d=d[d.geom_type=='Polygon'].copy();d['gcj_part']=d.groupby('block_id').cumcount()+1;d.insert(0,'feature_id',np.arange(1,len(d)+1,dtype='int32'));d=d.rename(columns={'area_km2':'source_area_km2'});d.crs=None
with warnings.catch_warnings():
 warnings.filterwarnings('ignore',message="'crs' was not provided.*");d.to_file(p,driver='GeoJSON')
print('features',len(d),'invalid',int((~d.is_valid).sum()),flush=True)
