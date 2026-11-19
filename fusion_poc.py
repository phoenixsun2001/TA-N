from pathlib import Path
import json, heapq
import numpy as np, geopandas as gpd, rasterio
from affine import Affine
from rasterio.warp import reproject,Resampling
from rasterio.features import rasterize,shapes
from scipy.ndimage import gaussian_filter, maximum_filter
from shapely.geometry import shape,Point
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

out=Path('outputs/imagery_poc'); out.mkdir(parents=True,exist_ok=True)
poi=gpd.read_file('outputs/commercial_poi.geojson').to_crs(32651)
roads=gpd.read_file('data/raw/osm/gis_osm_roads_free_1.shp').to_crs(32651)
rail=gpd.read_file('data/raw/osm/gis_osm_railways_free_1.shp').to_crs(32651)
water=gpd.read_file('data/raw/osm/gis_osm_water_a_free_1.shp').to_crs(32651)
roads=roads[roads.fclass.isin(['motorway','motorway_link','trunk','trunk_link']) & (roads.bridge!='T') & (roads.tunnel!='T') & (roads.layer.fillna(0)==0)]
rail=rail[rail.fclass.isin(['rail']) & (rail.bridge!='T') & (rail.tunnel!='T') & (rail.layer.fillna(0)==0)]

def region(seeds,blocked,cost,max_cost=1600):
 label=np.zeros((40,40),dtype='int16'); dist=np.full((40,40),np.inf); h=[]
 for i,(r,c) in enumerate(seeds,1):
  if not blocked[r,c]: label[r,c]=i; dist[r,c]=0; heapq.heappush(h,(0.,i,r,c))
 while h:
  d,i,r,c=heapq.heappop(h)
  if d>dist[r,c]+1e-5 or label[r,c]!=i: continue
  for rr,cc in ((r-1,c),(r+1,c),(r,c-1),(r,c+1)):
   if not (0<=rr<40 and 0<=cc<40) or blocked[rr,cc]: continue
   nd=d+50*(cost[r,c]+cost[rr,cc])
   if nd<=max_cost and nd<dist[rr,cc]: dist[rr,cc]=nd;label[rr,cc]=i;heapq.heappush(h,(nd,i,rr,cc))
 return label

def polygrid(arr,tr,name):
 rows=[{'id':int(v),'geometry':shape(g)} for g,v in shapes(arr.astype('int16'),transform=tr)]
 df=gpd.GeoDataFrame(rows,crs=32651).dissolve('id',as_index=False);df['method']=name;return df
fig,axes=plt.subplots(3,4,figsize=(16,12)); summaries={}
for row,site in enumerate(['xujiahui','lujiazui','jinqiao']):
 a=np.load(f'imagery/{site}_imagery.npz'); tr=Affine(*a['transform'][:6]); b=rasterio.transform.array_bounds(40,40,tr)
 with rasterio.open('imagery/viirs_2023_native.tif') as src:
  night=np.full((40,40),np.nan,dtype='float32'); reproject(rasterio.band(src,1),night,src_transform=src.transform,src_crs=src.crs,src_nodata=src.nodata,dst_transform=tr,dst_crs='EPSG:32651',dst_nodata=np.nan,resampling=Resampling.nearest)
 valid_s2=~np.isin(a['SCL'],[0,1,3,8,9,10,11]) & np.isfinite(a['SCL'])
 # Baseline 05.10 product metadata: reflectance=(DN-1000)/10000.
 green=(a['B03']-1000)/10000;red=(a['B04']-1000)/10000;nir=(a['B08']-1000)/10000;swir=(a['B11']-1000)/10000
 def nd(x,y): return np.divide(x-y,x+y,out=np.full_like(x,np.nan),where=(x+y)>0)
 ndvi=nd(nir,red);mndwi=nd(green,swir);ndbi=nd(swir,nir)
 vv=10*np.log10(np.maximum(a['vv'],1e-10));vh=10*np.log10(np.maximum(a['vh'],1e-10));sarvalid=np.isfinite(vv)&np.isfinite(vh)
 # Optical water and SAR low backscatter jointly indicate water. Missing imagery does not imply water.
 water_image=valid_s2 & sarvalid & (mndwi>0.10) & (vv < -13)
 layers=list(roads.cx[b[0]:b[2],b[1]:b[3]].geometry)+list(rail.cx[b[0]:b[2],b[1]:b[3]].geometry)+list(water.cx[b[0]:b[2],b[1]:b[3]].geometry)
 blocked=rasterize(((g,1) for g in layers if g is not None and not g.is_empty),out_shape=(40,40),transform=tr,all_touched=True).astype(bool)
 pts=poi.cx[b[0]:b[2],b[1]:b[3]]; cnt=np.zeros((40,40),dtype='float32')
 for p in pts.geometry:
  c=int((p.x-tr.c)/100);r=int((tr.f-p.y)/100)
  if 0<=r<40 and 0<=c<40:cnt[r,c]+=1
 den=gaussian_filter(cnt,3);cand=(den==maximum_filter(den,size=9)) & (den>=.08) & ~blocked
 seeds=[]
 for r,c in sorted(zip(*np.where(cand)),key=lambda p:-den[p]):
  if all((r-y)**2+(c-x)**2>=64 for y,x in seeds):seeds.append((int(r),int(c)))
 base=region(seeds,blocked,np.ones((40,40)))
 finite=night[np.isfinite(night)]; lo,hi=np.percentile(finite,[30,70]); night_cost=np.where(np.isfinite(night)&(night<lo),1.15,np.where(np.isfinite(night)&(night>hi),0.90,1.0))
 optical_cost=np.where(valid_s2&(ndvi>.35),1.8,1)
 sar_optical_cost=np.where(sarvalid & (vv>-8) & valid_s2 & (ndbi>0),.90,1)
 cost=night_cost*optical_cost*sar_optical_cost
 night_only=region(seeds,blocked,night_cost)
 optical_sar_only=region(seeds,blocked|water_image,optical_cost*sar_optical_cost)
 fused=region(seeds,blocked|water_image,cost)
 changed=(base!=fused);features=np.stack([night,ndvi,mndwi,ndbi,vv,vh,valid_s2.astype('float32'),water_image.astype('float32'),den,cost,base,fused])
 with rasterio.open(out/f'{site}_features.tif','w',driver='GTiff',height=40,width=40,count=len(features),dtype='float32',crs='EPSG:32651',transform=tr,nodata=-9999,compress='deflate') as ds:
  ds.write(np.where(np.isfinite(features),features,-9999).astype('float32'))
  for i,n in enumerate(['viirs_index','ndvi','mndwi','ndbi','s1_vv_db','s1_vh_db','s2_valid','image_water_candidate','poi_density','fusion_cost','baseline_id','fused_id'],1):ds.set_band_description(i,n)
 for label,name in [(base,'baseline'),(fused,'fused')]: polygrid(label,tr,name).to_crs(4326).to_file(out/f'{site}_{name}.geojson',driver='GeoJSON')
 gpd.GeoDataFrame([{'id':i,'geometry':Point(tr.c+(c+.5)*100,tr.f-(r+.5)*100)} for i,(r,c) in enumerate(seeds,1)],crs=32651).to_crs(4326).to_file(out/f'{site}_centers.geojson',driver='GeoJSON')
 summary={'poi_count':len(pts),'centers':len(seeds),'s2_valid_pct':float(valid_s2.mean()*100),'s1_valid_pct':float(sarvalid.mean()*100),'viirs_valid_pct':float(np.isfinite(night).mean()*100),'night_index_p30':float(lo),'night_index_p70':float(hi),'image_water_cells':int(water_image.sum()),'baseline_assigned_cells':int((base>0).sum()),'fused_assigned_cells':int((fused>0).sum()),'changed_cells':int(changed.sum()),'changed_km2':float(changed.sum()*.01),'night_only_changed_cells':int((night_only!=base).sum()),'s1_s2_only_changed_cells':int((optical_sar_only!=base).sum())};summaries[site]=summary
 for ax,im,title in [(axes[row,0],night,'VIIRS 2023 index'),(axes[row,1],mndwi,'S2 MNDWI'),(axes[row,2],base,'POI + OSM baseline'),(axes[row,3],fused,'With VIIRS + S1 + S2')]:
  ax.imshow(im,cmap='magma' if title.startswith('VIIRS') else 'RdBu' if title.startswith('S2') else 'tab20',interpolation='nearest');ax.set_title(f'{site} | {title}',fontsize=10);ax.set_xticks([]);ax.set_yticks([])
fig.tight_layout();fig.savefig(out/'comparison.png',dpi=160);plt.close(fig)
(out/'metrics.json').write_text(json.dumps(summaries,indent=2,ensure_ascii=False));print(json.dumps(summaries,indent=2))
