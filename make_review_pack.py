from pathlib import Path
import csv,json
import geopandas as gpd,numpy as np,rasterio
from rasterio.transform import rowcol,xy
from pyproj import Transformer
out=Path('outputs/city_imagery');rng=np.random.default_rng(20261114)
sites={'xujiahui':(121.4365,31.195),'lujiazui':(121.505,31.237),'jinqiao':(121.603,31.252)}
proj=Transformer.from_crs(4326,32651,always_xy=True);rev=Transformer.from_crs(32651,4326,always_xy=True)
rows=[]
with rasterio.open(out/'city_features.tif') as d:
 bands={name:d.read(i+1) for i,name in enumerate(d.descriptions)};base=bands['baseline_id'];fused=bands['fused_id'];water=bands['imagery_water_candidate'];s2=bands['s2_valid_dates'];s1=bands['s1_valid_dates'];ntl=bands['viirs_index']
 for site,lonlat in sites.items():
  x,y=proj.transform(*lonlat);r,c=rowcol(d.transform,x,y);window=np.zeros(base.shape,dtype=bool);window[max(0,r-8):min(d.height,r+8),max(0,c-8):min(d.width,c+8)]=True
  strata={'changed':window&(base!=fused)&(base>=0)&(fused>=0),'unchanged_assigned':window&(base==fused)&(base>0),'image_water':window&(water==1),'unassigned':window&(fused==0)&(water==0)}
  for stratum,condition in strata.items():
   indices=np.argwhere(condition); rng.shuffle(indices)
   for rr,cc in indices[:5]:
    xx,yy=xy(d.transform,int(rr),int(cc));lon,lat=rev.transform(xx,yy)
    rows.append({'sample_id':f'{site}_{stratum}_{len(rows)+1:03d}','site':site,'sample_type':stratum,'longitude_wgs84':round(lon,7),'latitude_wgs84':round(lat,7),'baseline_id':int(base[rr,cc]),'fused_id':int(fused[rr,cc]),'s2_valid_dates':int(s2[rr,cc]),'s1_valid_dates':int(s1[rr,cc]),'viirs_index':float(ntl[rr,cc]),'manual_label':'','evidence':'','notes':''})
with open(out/'review_samples.csv','w',newline='',encoding='utf-8-sig') as f:
 writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
print('review samples',len(rows))
frame=gpd.GeoDataFrame(rows,geometry=gpd.points_from_xy([r['longitude_wgs84'] for r in rows],[r['latitude_wgs84'] for r in rows]),crs=4326)
frame.to_file(out/'review_samples.geojson',driver='GeoJSON')
