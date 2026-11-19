from pathlib import Path
import json,warnings
import numpy as np,geopandas as g,shapely,rasterio
from rasterio.features import rasterize
from shapely.geometry import Polygon,Point
from shapely.strtree import STRtree
from convert_gcj02 import convert
out=Path('outputs/improve_v7');out.mkdir(exist_ok=True)
d=g.read_file('outputs/rebuild_v6/reconstructed_blocks.gpkg').to_crs(32651);geom=d.geometry.to_numpy();tree=STRtree(geom)
water=g.read_file('data/raw/osm/gis_osm_water_a_free_1.shp').to_crs(32651);water=water[water.fclass.isin(['water','riverbank','reservoir','dock'])];wg=shapely.make_valid(water.geometry.array);wt=STRtree(wg)
import pandas as pd
roads=g.read_file('data/raw/osm/gis_osm_roads_free_1.shp').to_crs(32651);extra=g.read_file('data/user_shanghai/national_roads_candidate_filtered.gpkg').to_crs(32651);roads=g.GeoDataFrame(pd.concat([roads,extra],ignore_index=True),crs=32651)
width={'motorway':15,'motorway_link':8,'trunk':12,'trunk_link':7,'primary':9,'primary_link':6,'secondary':7,'secondary_link':5,'tertiary':6,'tertiary_link':5}
roads=roads[roads.fclass.isin(width)&(roads.tunnel!='T')&(roads.layer.fillna(0)>=0)];rb=shapely.buffer(roads.geometry.array,roads.fclass.map(width).to_numpy(),quad_segs=2);rt=STRtree(rb)
rail=g.read_file('data/raw/osm/gis_osm_railways_free_1.shp').to_crs(32651);rail=rail[(rail.fclass.isin(['rail','light_rail','tram','monorail','narrow_gauge','funicular'])|((rail.fclass=='subway')&((rail.bridge=='T')|(rail.layer.fillna(0)>0))))&(rail.tunnel!='T')&(rail.layer.fillna(0)>=0)];rr=shapely.buffer(rail.geometry.array,5,quad_segs=2);rlt=STRtree(rr)
# Small holes are filled only where no mapped water or other block occupies them.
filled=protected=0;geoms=[]
for i,p in enumerate(geom):
 holes=[]
 for ring in p.interiors:
  h=Polygon(ring)
  if h.area>1000:holes.append(ring);continue
  near=tree.query(h,predicate='intersects');others=near[near!=i]
  occupied=any(h.intersection(geom[j]).area>.01 for j in others)
  water_ix=wt.query(h,predicate='intersects');water_overlap=any(h.intersection(wg[j]).area>.01 for j in water_ix)
  road_overlap=any(h.intersection(rb[j]).area>.01 for j in rt.query(h,predicate='intersects'))
  rail_overlap=any(h.intersection(rr[j]).area>.01 for j in rlt.query(h,predicate='intersects'))
  if occupied or water_overlap or road_overlap or rail_overlap:holes.append(ring);protected+=1
  else:filled+=1
 geoms.append(Polygon(p.exterior,holes))
d.geometry=geoms;d['area_km2']=d.area/1e6;d['holes_remaining']=[len(p.interiors) for p in d.geometry]
with rasterio.open('outputs/city_imagery/city_features.tif') as ds:
 cell=rasterize(((z,i+1) for i,z in enumerate(d.geometry)),out_shape=(ds.height,ds.width),transform=ds.transform,fill=0,dtype='int32')-1
 for name in ['viirs_index','ndvi','ndbi','s1_vv_db']:
  v=ds.read(ds.descriptions.index(name)+1);ok=(cell>=0)&(v!=-9999)&np.isfinite(v);counts=np.bincount(cell[ok],minlength=len(d));sums=np.bincount(cell[ok],weights=v[ok],minlength=len(d));d[name+'_mean']=np.divide(sums,counts,out=np.full(len(d),np.nan),where=counts>0)
 d['imagery_pixels']=np.bincount(cell[cell>=0],minlength=len(d))
land=g.read_file('data/raw/osm/gis_osm_landuse_a_free_1.shp').to_crs(32651);land=land[land.fclass.isin(['residential','industrial','farmland','forest','park','commercial','retail','grass','cemetery'])].copy()
centers=g.GeoDataFrame({'block_id':d.block_id},geometry=d.geometry.representative_point(),crs=32651);join=g.sjoin(centers,land[['fclass','geometry']],predicate='within').drop_duplicates('block_id');d['landuse_center']=d.block_id.map(join.set_index('block_id').fclass)
d['review_fragment']=(d.area<1000);d['review_large']=(d.area>2e6);d['compactness']=np.divide(4*np.pi*d.area,d.length**2,out=np.zeros(len(d)),where=d.length>0)
d.to_file(out/'blocks_projected.gpkg',driver='GPKG');gcj=d.to_crs(4326);gcj.geometry=shapely.transform(gcj.geometry.array,convert,include_z=False);bad=~gcj.is_valid;gcj.loc[bad,'geometry']=shapely.make_valid(gcj.loc[bad,'geometry'].array);gcj=gcj.explode(index_parts=False,ignore_index=True);gcj=gcj[gcj.geom_type=='Polygon'];gcj.crs=None
with warnings.catch_warnings():
 warnings.filterwarnings('ignore',message="'crs' was not provided.*");gcj.to_file(out/'blocks_gcj02.geojson',driver='GeoJSON')
m={'blocks':len(d),'holes_filled_le_1000m2':filled,'small_holes_protected_water_or_blocks':protected,'remaining_holes':int(d.holes_remaining.sum()),'invalid':int((~d.is_valid).sum()),'imagery_zero_pixel_blocks':int((d.imagery_pixels==0).sum()),'gcj_features':len(gcj),'gcj_invalid':int((~gcj.is_valid).sum())};(out/'metrics.json').write_text(json.dumps(m,indent=2));print(json.dumps(m,indent=2))
