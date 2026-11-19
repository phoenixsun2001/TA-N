from pathlib import Path
import json, heapq, hashlib
import numpy as np, geopandas as gpd, rasterio
from scipy.ndimage import gaussian_filter, maximum_filter
from rasterio.features import rasterize, shapes
from rasterio.transform import from_origin
from shapely.geometry import shape, Point
from shapely import make_valid
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parent; raw=ROOT/'data/raw'; out=ROOT/'outputs'; out.mkdir(exist_ok=True)
CRS='EPSG:32651'; RES=250; MAX_DIST=2000
classes={'restaurant','cafe','fast_food','bar','pub','bakery','food_court','ice_cream','convenience','supermarket','clothes','shoe_shop','department_store','mall','marketplace','greengrocer','beverages','jeweller','bookshop','beauty_shop','hairdresser','cinema','theatre','nightclub','fitness_centre','sports_centre','gift_shop','general','kiosk','toy_shop','furniture_shop','mobile_phone_shop','pharmacy'}
def read(name): return gpd.read_file(raw/'osm'/f'gis_osm_{name}_free_1.shp').to_crs(CRS)
district=gpd.read_file(raw/'shanghai_boundary.geojson').to_crs(CRS); city=make_valid(district.geometry.union_all())
minx,miny,maxx,maxy=city.bounds; minx=np.floor(minx/RES)*RES; maxy=np.ceil(maxy/RES)*RES
width=int(np.ceil((maxx-minx)/RES)); height=int(np.ceil((maxy-miny)/RES)); transform=from_origin(minx,maxy,RES,RES)
mask=rasterize([(city,1)],out_shape=(height,width),transform=transform,all_touched=True).astype(bool)
p=read('pois'); p=p[p.fclass.isin(classes) & p.geometry.within(city)].copy()
xy=np.array([(v.x,v.y) for v in p.geometry]); cols=((xy[:,0]-minx)/RES).astype(int); rows=((maxy-xy[:,1])/RES).astype(int)
counts=np.zeros(mask.shape,dtype='float32'); np.add.at(counts,(rows,cols),1)
density=gaussian_filter(counts, sigma=1.2)
water=read('water_a'); water=water[water.intersects(city)]
watermask=rasterize(((x,1) for x in water.geometry if x is not None and not x.is_empty),out_shape=mask.shape,transform=transform,all_touched=False).astype(bool)
roads=read('roads'); roads=roads[roads.fclass.isin(['motorway','motorway_link','trunk','trunk_link']) & (roads.bridge!='T') & (roads.tunnel!='T') & (roads.layer.fillna(0)==0) & roads.intersects(city)]
rail=read('railways'); rail=rail[rail.fclass.isin(['rail']) & (rail.bridge!='T') & (rail.tunnel!='T') & (rail.layer.fillna(0)==0) & rail.intersects(city)]
# Coarse raster barriers are deliberately approximate; grid cells intersecting lines are blocked.
barrier=rasterize(((x,1) for x in list(roads.geometry)+list(rail.geometry) if x is not None and not x.is_empty),out_shape=mask.shape,transform=transform,all_touched=True).astype(bool)
blocked=watermask|barrier
candidate=(density==maximum_filter(density,size=5)) & (density>=0.40) & mask & ~blocked
indices=sorted(zip(*np.where(candidate)),key=lambda rc:-density[rc]); seeds=[]
for r,c in indices:
 if all((r-a)**2+(c-b)**2>=16 for a,b in seeds): seeds.append((r,c))
labels=np.zeros(mask.shape,dtype='int32'); dist=np.full(mask.shape,np.inf,dtype='float32'); heap=[]
for i,(r,c) in enumerate(seeds,1): labels[r,c]=i; dist[r,c]=0; heapq.heappush(heap,(0,i,int(r),int(c)))
# 4-neighbour region growing; distances are grid cost, not walking times.
while heap:
 d,i,r,c=heapq.heappop(heap)
 if d>dist[r,c]+0.01 or labels[r,c]!=i: continue
 for rr,cc in ((r-1,c),(r+1,c),(r,c-1),(r,c+1)):
  nd=d+RES
  if 0<=rr<height and 0<=cc<width and mask[rr,cc] and not blocked[rr,cc] and nd<=MAX_DIST and nd<dist[rr,cc]:
   dist[rr,cc]=nd; labels[rr,cc]=i; heapq.heappush(heap,(nd,i,rr,cc))
# Each cell has exactly one label. Negative values represent non-assigned / proxy barriers.
zone=labels.copy(); zone[mask & (zone==0)]=-1; zone[mask & barrier & ~watermask]=-2; zone[mask & watermask]=-3; zone[~mask]=0
core=np.where((density>=0.40)&(labels>0)&~blocked,labels,0).astype('int32')
def polygons(arr, include):
 features=[]
 for geom,val in shapes(arr,mask=include,transform=transform):
  geom=make_valid(shape(geom).intersection(city))
  if not geom.is_empty and geom.area>0: features.append({'id':int(val),'geometry':geom})
 return gpd.GeoDataFrame(features,crs=CRS).dissolve(by='id',as_index=False)
z=polygons(zone,zone!=0); z['zone_type']=z.id.map(lambda x: 'commercial_candidate' if x>0 else {-1:'unassigned',-2:'barrier_proxy',-3:'water_proxy'}[x]); z['method']='250m_grid_cost_not_walking'; z['area_km2']=z.area/1e6
cores=polygons(core,core>0); cores['area_km2']=cores.area/1e6
centers=gpd.GeoDataFrame([{'id':i,'density':float(density[r,c]),'geometry':Point(minx+(c+.5)*RES,maxy-(r+.5)*RES)} for i,(r,c) in enumerate(seeds,1)],crs=CRS)
for df,name in [(z,'trade_area_partition'),(cores,'commercial_core'),(centers,'centers'),(p,'commercial_poi')]: df.to_crs('EPSG:4326').to_file(out/f'{name}.geojson',driver='GeoJSON')
# The city boundary and the original partition are in the same projection; calculate topology after clipping.
union=z.geometry.union_all(); coverage_gap=city.difference(union).area; overlap=z.area.sum()-union.area
info={'grid_m':RES,'max_grid_distance_m':MAX_DIST,'density_sigma_cells':1.2,'density_threshold':0.40,'center_count':len(centers),'poi_count':len(p),'core_count':len(cores),'city_area_km2':city.area/1e6,'gap_m2':coverage_gap,'overlap_m2':overlap,'invalid_geometries':int((~z.is_valid).sum()),'zone_area_km2':z.groupby('zone_type').area_km2.sum().to_dict(),'input_sha256':{str(f.relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest() for f in [raw/'shanghai.zip',raw/'shanghai_boundary.geojson']}}
(out/'metrics.json').write_text(json.dumps(info,ensure_ascii=False,indent=2))
fig,ax=plt.subplots(figsize=(12,10)); z[z.id>0].plot(ax=ax,column='id',cmap='tab20',linewidth=0,alpha=.7); z[z.id==-1].plot(ax=ax,color='#e5e7eb',linewidth=0); z[z.id==-2].plot(ax=ax,color='#8f959f',linewidth=0); z[z.id==-3].plot(ax=ax,color='#a6cee3',linewidth=0); cores.boundary.plot(ax=ax,color='#e44b29',linewidth=.35); ax.set_axis_off(); ax.set_title('Shanghai commercial area PoC | OSM | 250 m grid',fontsize=15); fig.tight_layout(); fig.savefig(out/'overview.png',dpi=180); plt.close(fig)
print(json.dumps(info,ensure_ascii=False,indent=2))
