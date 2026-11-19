from pathlib import Path
import numpy as np,geopandas as g,shapely,json
from shapely.geometry import box
out=Path('outputs/revised_poc');out.mkdir(exist_ok=True)
b=g.read_file('outputs/block_fit/city_blocks.gpkg').to_crs(32651);b['old_id']=b.id;b['synthetic_split']=False;b['active']=~((~b.is_water)&(b.area>2e6));b['geometry']=shapely.make_valid(b.geometry.array)
large=b[~b.active];print('oversized land blocks',len(large),flush=True)
new=[];step=1000
for row in large.itertuples():
 minx,miny,maxx,maxy=row.geometry.bounds
 xs=np.arange(np.floor(minx/step)*step,maxx,step);ys=np.arange(np.floor(miny/step)*step,maxy,step)
 X,Y=np.meshgrid(xs,ys);grid=shapely.box(X.ravel(),Y.ravel(),X.ravel()+step,Y.ravel()+step);hit=shapely.intersects(row.geometry,grid);pieces=shapely.intersection(row.geometry,grid[hit])
 for piece in pieces:
  for poly in shapely.get_parts(piece):
   if poly.geom_type=='Polygon' and poly.area>.00001:new.append({'old_id':int(row.id),'is_water':False,'synthetic_split':True,'active':True,'geometry':poly})
print('fallback units',len(new),flush=True)
import pandas as pd
children=g.GeoDataFrame(new,crs=32651);b=g.GeoDataFrame(pd.concat([b,children],ignore_index=True),crs=32651);b['id']=np.arange(len(b),dtype='int32');print('active',int(b.active.sum()),flush=True);b.to_file(out/'activity_blocks.gpkg',driver='GPKG');(out/'refinement_metrics.json').write_text(json.dumps({'original_blocks':len(b)-len(children),'large_land_blocks_split':len(large),'fallback_pieces':len(children),'active_faces':int(b.active.sum())},indent=2));print('done')
