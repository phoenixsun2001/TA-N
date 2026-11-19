from pathlib import Path
import json
import numpy as np
import rasterio
from rasterio.warp import transform_bounds, reproject, Resampling
from rasterio.transform import from_origin
from rasterio.windows import from_bounds, transform as window_transform
from pyproj import Transformer
from pystac_client import Client
import planetary_computer as pc

OUT=Path('imagery'); OUT.mkdir(exist_ok=True)
# 4 x 4 km pilot squares, centers are approximate and not administrative boundaries.
sites={'xujiahui':(121.4365,31.195),'lujiazui':(121.505,31.237),'jinqiao':(121.603,31.252)}
catalog=Client.open('https://planetarycomputer.microsoft.com/api/stac/v1')
s2=catalog.get_collection('sentinel-2-l2a').get_item('S2A_MSIL2A_20231224T024121_R089_T51RUQ_20231224T075251')
s1ids=['S1A_IW_GRDH_1SDV_20231229T095512_20231229T095537_051868_064426_rtc','S1A_IW_GRDH_1SDV_20231217T095513_20231217T095538_051693_063E16_rtc','S1A_IW_GRDH_1SDV_20231205T095513_20231205T095538_051518_0637F8_rtc']
s1=[catalog.get_collection('sentinel-1-rtc').get_item(i) for i in s1ids]
project=Transformer.from_crs(4326,32651,always_xy=True)
def fetch(url,tr,size,resampling):
 with rasterio.Env(GDAL_DISABLE_READDIR_ON_OPEN='EMPTY_DIR',GDAL_HTTP_TIMEOUT='60',GDAL_HTTP_MAX_RETRY='3'):
  with rasterio.open(url) as src:
   bounds=rasterio.transform.array_bounds(size,size,tr); bounds=transform_bounds('EPSG:32651',src.crs,*bounds,densify_pts=21)
   window=from_bounds(*bounds,transform=src.transform).round_offsets().round_lengths(); data=src.read(1,window=window,masked=True)
   dest=np.full((size,size),np.nan,dtype='float32')
   reproject(data.astype('float32').filled(np.nan),dest,src_transform=window_transform(window,src.transform),src_crs=src.crs,src_nodata=np.nan,dst_transform=tr,dst_crs='EPSG:32651',dst_nodata=np.nan,resampling=resampling)
   return dest
for name,xy in sites.items():
 target=OUT/f'{name}_imagery.npz'
 if target.exists(): print(name,'cached',flush=True); continue
 x,y=project.transform(*xy); x=round(x/100)*100; y=round(y/100)*100; tr=from_origin(x-2000,y+2000,100,100); size=40
 result={}
 for band in ['B03','B04','B08','B11','SCL']:
  print(name,band,flush=True); result[band]=fetch(pc.sign(s2.assets[band].href),tr,size,Resampling.nearest if band=='SCL' else Resampling.average)
 for band in ['vv','vh']:
  layers=[]
  for item in s1:
   print(name,band,item.datetime,flush=True); layers.append(fetch(pc.sign(item.assets[band].href),tr,size,Resampling.average))
  result[band]=np.nanmedian(layers,axis=0)
 result['transform']=np.array(tuple(tr)); np.savez_compressed(target,**result)
meta={'s2_id':s2.id,'s2_datetime':str(s2.datetime),'s2_processing_baseline':s2.properties.get('s2:processing_baseline'),'s2_cloud_percent_scene':s2.properties.get('eo:cloud_cover'),'s1_ids':s1ids,'sites':sites,'grid_m':100,'site_size_km':4}
(OUT/'sources.json').write_text(json.dumps(meta,ensure_ascii=False,indent=2))
