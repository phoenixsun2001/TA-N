from pathlib import Path
import json, heapq, warnings
import numpy as np, geopandas as gpd, rasterio
from affine import Affine
from rasterio.features import rasterize,shapes
from rasterio.warp import reproject,Resampling
from scipy.ndimage import gaussian_filter
from shapely.geometry import shape
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
warnings.filterwarnings('ignore',message='All-NaN slice encountered')
cache=Path('imagery/city_cache');out=Path('outputs/city_imagery');out.mkdir(parents=True,exist_ok=True)
grid=np.load(cache/'grid.npz'); tr=Affine(*grid['transform'][:6]);mask=grid['mask'];h,w=mask.shape
meta=json.loads((cache/'sources.json').read_text());offsets=json.loads((cache/'s2_offsets.json').read_text()); assert {x['id'] for x in offsets}=={x['id'] for x in meta['sentinel2']} and all(x['offsets']==['-1000'] and x['quantification']==['10000'] for x in offsets), 'S2 offsets missing or differ; verify scene metadata before processing';boundary=gpd.read_file('data/raw/shanghai_boundary.geojson').to_crs(32651).geometry.union_all()
def load(id,band):return np.load(cache/f'{id}_{band}.npy')
def norm(x,y):return np.divide(x-y,x+y,out=np.full_like(x,np.nan),where=np.isfinite(x)&np.isfinite(y)&(x>=0)&(y>=0)&((x+y)>1e-5))
# First merge overlapping same-date tiles; then median across independent dates, preserving each pixel's valid count.
s2dates=sorted({a['datetime'][:10] for a in meta['sentinel2']});daily=[]
for date in s2dates:
 maps=[np.full((h,w),np.nan,dtype='float32') for _ in range(3)]
 for info in [x for x in meta['sentinel2'] if x['datetime'].startswith(date)]:
  id=info['id'];scl=load(id,'SCL');valid=np.isfinite(scl)&~np.isin(scl,[0,1,3,7,8,9,10,11]); bands=[(load(id,b)-1000)/10000 for b in ['B03','B04','B08','B11']];green,red,nir,swir=bands
  for dest,vals in zip(maps,[norm(nir,red),norm(green,swir),norm(swir,nir)]):
   use=valid&np.isfinite(vals)&~np.isfinite(dest);dest[use]=vals[use]
 daily.append(maps)
stack=np.array(daily);ndvi,mndwi,ndbi=np.nanmedian(stack,axis=0);s2obs=np.isfinite(stack[:,0]).sum(axis=0)
# Merge frames from each same-orbit overpass before taking the temporal median.
s1dates=sorted({a['datetime'][:10] for a in meta['sentinel1']});sardaily=[]
for date in s1dates:
 maps=[np.full((h,w),np.nan,dtype='float32') for _ in range(2)]
 for info in [x for x in meta['sentinel1'] if x['datetime'].startswith(date)]:
  for dest,band in zip(maps,['vv','vh']):
   v=load(info['id'],band);use=np.isfinite(v)&(v>0)&~np.isfinite(dest);dest[use]=v[use]
 sardaily.append(maps)
sarstack=np.array(sardaily);sarobs=np.isfinite(sarstack[:,0]).sum(axis=0);vvp,vhp=np.nanmedian(sarstack,axis=0);vv=10*np.log10(np.maximum(vvp,1e-10));vh=10*np.log10(np.maximum(vhp,1e-10))
with rasterio.open(cache/'viirs_2023_native.tif') as ds:
 night=np.full((h,w),np.nan,dtype='float32');reproject(rasterio.band(ds,1),night,src_transform=ds.transform,src_crs=ds.crs,src_nodata=ds.nodata,dst_transform=tr,dst_crs='EPSG:32651',dst_nodata=np.nan,resampling=Resampling.nearest)
night[~mask]=np.nan
s2valid=s2obs>0;sarvalid=sarobs>0;watimg=mask&s2valid&sarvalid&(mndwi>.10)&(vv < -13)
def read_osm(name):return gpd.read_file(f'data/raw/osm/gis_osm_{name}_free_1.shp').to_crs(32651)
roads=read_osm('roads');roads=roads[roads.fclass.isin(['motorway','motorway_link','trunk','trunk_link'])&(roads.bridge!='T')&(roads.tunnel!='T')&(roads.layer.fillna(0)==0)&roads.intersects(boundary)]
rail=read_osm('railways');rail=rail[rail.fclass.isin(['rail'])&(rail.bridge!='T')&(rail.tunnel!='T')&(rail.layer.fillna(0)==0)&rail.intersects(boundary)]
water=read_osm('water_a');water=water[water.intersects(boundary)]
barrier=rasterize(((g,1) for g in list(roads.geometry)+list(rail.geometry)),out_shape=(h,w),transform=tr,all_touched=True).astype(bool)
watermask=rasterize(((g,1) for g in water.geometry),out_shape=(h,w),transform=tr,all_touched=False).astype(bool);blocked=barrier|watermask|~mask
centers=gpd.read_file('outputs/centers.geojson').to_crs(32651);seeds=[]
for p in centers.itertuples():
 r=int((tr.f-p.geometry.y)/250);c=int((p.geometry.x-tr.c)/250)
 if 0<=r<h and 0<=c<w and not blocked[r,c]:seeds.append((int(p.id),r,c))
print('seeds',len(seeds),flush=True)
def region(block,cost):
 label=np.zeros((h,w),dtype='int16');dist=np.full((h,w),np.inf,dtype='float32');q=[]
 for i,r,c in seeds:
  if not block[r,c]:label[r,c]=i;dist[r,c]=0;heapq.heappush(q,(0.,i,r,c))
 while q:
  d,i,r,c=heapq.heappop(q)
  if d>dist[r,c]+.01 or label[r,c]!=i:continue
  for rr,cc in ((r-1,c),(r+1,c),(r,c-1),(r,c+1)):
   if not (0<=rr<h and 0<=cc<w) or block[rr,cc]:continue
   nd=d+125*(float(cost[r,c])+float(cost[rr,cc]))
   if nd<=2000 and nd<dist[rr,cc]-0.001:dist[rr,cc]=nd;label[rr,cc]=i;heapq.heappush(q,(nd,i,rr,cc))
 return label
lights=night[mask&np.isfinite(night)];lo,hi=np.percentile(lights,[30,70]);nightcost=np.where(np.isfinite(night)&(night<lo),1.15,np.where(np.isfinite(night)&(night>hi),.90,1.0));optcost=np.where(s2valid&(ndvi>.35),1.8,1.0);saropt=np.where(sarvalid&s2valid&(vv>-8)&(ndbi>0),.90,1.0);cost=(nightcost*optcost*saropt).astype('float32')
base=region(blocked,np.ones((h,w),dtype='float32'));night_only=region(blocked,nightcost);sentinel_only=region(blocked|watimg,optcost*saropt);fused=region(blocked|watimg,cost)
def classified(a,include_image=False):
 z=a.copy();z[mask&(z==0)]=-1;z[mask&barrier&~watermask]=-2;z[mask&watermask]=-3
 if include_image:z[mask&watimg&~watermask&~barrier]=-4
 z[~mask]=0;return z
bz=classified(base);fz=classified(fused,True)
poi=gpd.read_file('outputs/commercial_poi.geojson').to_crs(32651);counts=np.zeros((h,w),dtype='float32')
for p in poi.geometry:
 r=int((tr.f-p.y)/250);c=int((p.x-tr.c)/250)
 if 0<=r<h and 0<=c<w:counts[r,c]+=1
density=gaussian_filter(counts,1.2);core=np.where((density>=.4)&(fused>0)&~(blocked|watimg),fused,0).astype('int16')
features=np.stack([night,ndvi,mndwi,ndbi,vv,vh,s2obs,sarobs,watimg,density,cost,base,fused]);names=['viirs_index','ndvi','mndwi','ndbi','s1_vv_db','s1_vh_db','s2_valid_dates','s1_valid_dates','imagery_water_candidate','poi_density','fusion_cost','baseline_id','fused_id']
with rasterio.open(out/'city_features.tif','w',driver='GTiff',height=h,width=w,count=len(names),dtype='float32',crs=32651,transform=tr,nodata=-9999,compress='deflate') as ds:
 features[:,~mask]=np.nan;ds.write(np.where(np.isfinite(features),features,-9999).astype('float32'))
 for i,n in enumerate(names,1):ds.set_band_description(i,n)
def polygonize(a,path):
 rows=[]
 for geom,value in shapes(a.astype('int16'),mask=a!=0,transform=tr):
  part=shape(geom).intersection(boundary)
  if not part.is_empty and part.area>0:rows.append({'id':int(value),'geometry':part})
 df=gpd.GeoDataFrame(rows,crs=32651).dissolve(by='id',as_index=False);df['zone_type']=df.id.map(lambda x:'candidate' if x>0 else {-1:'unassigned',-2:'road_rail_proxy',-3:'osm_water_proxy',-4:'imagery_water_candidate'}.get(x,'unknown'));df['area_km2']=df.area/1e6;df.to_crs(4326).to_file(path,driver='GeoJSON');return df
print('polygonizing',flush=True)
b=polygonize(bz,out/'city_baseline.geojson');f=polygonize(fz,out/'city_fused.geojson');co=polygonize(core,out/'city_commercial_core.geojson')
union=f.geometry.union_all();total=mask.sum();pct=lambda a:float((a&mask).sum()/total*100)
metrics={'grid_m':250,'input_area_km2':float(boundary.area/1e6),'poi_count':len(poi),'seed_count':len(seeds),'s2_dates':s2dates,'s1_dates':s1dates,'s2_coverage_pct':pct(s2valid),'s2_two_dates_pct':pct(s2obs>=2),'s1_coverage_pct':pct(sarvalid),'s1_three_dates_pct':pct(sarobs>=3),'viirs_coverage_pct':pct(np.isfinite(night)),'viirs_index_p30':float(lo),'viirs_index_p70':float(hi),'image_water_cells':int(watimg.sum()),'changed_cells':int(((base!=fused)&mask).sum()),'night_only_changed_cells':int(((base!=night_only)&mask).sum()),'s1_s2_only_changed_cells':int(((base!=sentinel_only)&mask).sum()),'grid_changed_km2':float(((base!=fused)&mask).sum()*0.0625),'baseline_assigned_grid_pct':pct(base>0),'fused_assigned_grid_pct':pct(fused>0),'topology_gap_m2':float(boundary.difference(union).area),'topology_overlap_m2':float(f.area.sum()-union.area),'invalid_geometries':int((~f.is_valid).sum()),'fused_zone_areas_km2':f.groupby('zone_type').area_km2.sum().to_dict()}
(out/'metrics.json').write_text(json.dumps(metrics,ensure_ascii=False,indent=2));print(json.dumps(metrics,indent=2),flush=True)
fig,axs=plt.subplots(2,3,figsize=(18,12))
for ax,arr,title,cmap in [(axs[0,0],night,'VIIRS 2023 relative index','magma'),(axs[0,1],ndvi,'Sentinel-2 NDVI median','RdYlGn'),(axs[0,2],vv,'Sentinel-1 VV dB median','gray'),(axs[1,0],np.where(mask,bz,np.nan),'POI + OSM baseline','tab20'),(axs[1,1],np.where(mask,fz,np.nan),'Satellite fusion','tab20'),(axs[1,2],np.where(mask,(base!=fused).astype(float),np.nan),'Assignment changed','Reds')]:ax.imshow(arr,cmap=cmap,interpolation='nearest');ax.set_title(title);ax.axis('off')
fig.tight_layout();fig.savefig(out/'city_overview.png',dpi=160);plt.close(fig)
