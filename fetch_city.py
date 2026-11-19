from pathlib import Path
import json,re
import numpy as np, geopandas as gpd, rasterio
from rasterio.vrt import WarpedVRT
from rasterio.enums import Resampling
from rasterio.transform import from_origin
from rasterio.features import geometry_mask
from shapely.geometry import shape
from pyproj import Transformer
from pystac_client import Client
import planetary_computer as pc

OUT=Path('imagery/city_cache');OUT.mkdir(parents=True,exist_ok=True)
city=gpd.read_file('data/raw/shanghai_boundary.geojson').to_crs(32651).geometry.union_all();x0,y0,x1,y1=city.bounds; res=250;x0=np.floor(x0/res)*res;y1=np.ceil(y1/res)*res;w=int(np.ceil((x1-x0)/res));h=int(np.ceil((y1-y0)/res));tr=from_origin(x0,y1,res,res)
mask=geometry_mask([city],out_shape=(h,w),transform=tr,invert=True,all_touched=True)
np.savez_compressed(OUT/'grid.npz',transform=np.array(tuple(tr)),mask=mask)
print('grid',w,h, 'city cells',int(mask.sum()),flush=True)
c=Client.open('https://planetarycomputer.microsoft.com/api/stac/v1');bound=gpd.GeoSeries([city],crs=32651).to_crs(4326).iloc[0]; bbox=bound.bounds
items=list(c.search(collections=['sentinel-2-l2a'],bbox=bbox,datetime='2023-12-01/2023-12-25',max_items=500).items())
chosen=[]
for day in ['2023-12-24','2023-12-01']:
 tiles={}
 for item in items:
  if str(item.datetime.date())!=day or not shape(item.geometry).intersects(bound):continue
  tile=item.properties.get('s2:mgrs_tile'); cloud=item.properties.get('eo:cloud_cover',100)
  # prefer original processing among duplicates, then low cloud
  key=(cloud,item.id)
  if tile not in tiles or key<tiles[tile][0]:tiles[tile]=(key,item)
 chosen += [v[1] for v in tiles.values()]
print('S2 scenes',len(chosen),[x.id for x in chosen],flush=True)
s1s=list(c.search(collections=['sentinel-1-rtc'],bbox=bbox,datetime='2023-12-01/2023-12-31',max_items=100).items());s1s=[x for x in s1s if x.properties.get('sat:relative_orbit')==171 and shape(x.geometry).intersects(bound)]
print('S1 scenes',len(s1s),[x.id for x in s1s],flush=True)
meta={'crs':'EPSG:32651','grid_m':res,'shape':[h,w],'sentinel2':[{'id':x.id,'datetime':str(x.datetime),'cloud_scene_pct':x.properties.get('eo:cloud_cover'),'baseline':x.properties.get('s2:processing_baseline')} for x in chosen],'sentinel1':[{'id':x.id,'datetime':str(x.datetime),'orbit':x.properties.get('sat:relative_orbit')} for x in s1s]}
(OUT/'sources.json').write_text(json.dumps(meta,ensure_ascii=False,indent=2))
def fetch(item,band,method):
 target=OUT/f'{item.id}_{band}.npy'
 if target.exists():return
 print('fetch',item.id,band,flush=True)
 with rasterio.Env(GDAL_DISABLE_READDIR_ON_OPEN='EMPTY_DIR',GDAL_HTTP_TIMEOUT='90',GDAL_HTTP_MAX_RETRY='3',CPL_VSIL_CURL_CACHE_SIZE='20000000'):
  with rasterio.open(pc.sign(item.assets[band].href), **({'overview_level':3} if band in ['vv','vh'] else {})) as src:
   with WarpedVRT(src,crs='EPSG:32651',transform=tr,width=w,height=h,dtype='float32',nodata=-9999,resampling=method) as vrt:
    a=vrt.read(1,masked=True).filled(np.nan);a[~mask]=np.nan;np.save(target,a)
from concurrent.futures import ThreadPoolExecutor
tasks=[(i,b,Resampling.nearest if b=='SCL' else Resampling.average) for i in chosen for b in ['B03','B04','B08','B11','SCL']] + [(i,b,Resampling.average) for i in s1s for b in ['vv','vh']]
with ThreadPoolExecutor(max_workers=6) as pool:
 futures=[pool.submit(fetch,*task) for task in tasks]
 for future in futures: future.result()
url='https://zenodo.org/records/17057876/files/nightlights.average_viirs.v21_m_500m_s_20230101_20231231_go_epsg4326_v20250904.tif?download=1'
target=OUT/'viirs_2023_native.tif'
if not target.exists():
 from rasterio.windows import from_bounds
 with rasterio.Env(GDAL_DISABLE_READDIR_ON_OPEN='EMPTY_DIR',GDAL_HTTP_TIMEOUT='90'):
  with rasterio.open(url) as src:
   win=from_bounds(bbox[0]-.02,bbox[1]-.02,bbox[2]+.02,bbox[3]+.02,src.transform).round_offsets().round_lengths();a=src.read(1,window=win);p=src.profile.copy();p.update(width=int(win.width),height=int(win.height),transform=src.window_transform(win),compress='deflate')
   with rasterio.open(target,'w',**p) as dst:dst.write(a,1);dst.update_tags(source=url)
print('download complete',flush=True)
