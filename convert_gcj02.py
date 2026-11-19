from pathlib import Path
import json, warnings
import geopandas as gpd
import numpy as np, shapely

src=Path('outputs/block_fit');out=src/'gcj02';out.mkdir(exist_ok=True)
# Common open-source WGS84 -> GCJ-02 approximation; numpy vectorized.
def convert(c):
 if len(c)==0:return c
 lon=c[:,0];lat=c[:,1];x=lon-105.;y=lat-35.;p=np.pi
 d=20*np.sin(6*x*p)+20*np.sin(2*x*p)
 dlat=(d+20*np.sin(y*p)+40*np.sin(y*p/3)+160*np.sin(y*p/12)+320*np.sin(y*p/30))*2/3
 dlng=(d+20*np.sin(x*p)+40*np.sin(x*p/3)+150*np.sin(x*p/12)+300*np.sin(x*p/30))*2/3
 dlat+=-100+2*x+3*y+.2*y*y+.1*x*y+.2*np.sqrt(np.abs(x))
 dlng+=300+x+2*y+.1*x*x+.1*x*y+.1*np.sqrt(np.abs(x))
 rad=lat/180*p;magic=1-.00669342162296594323*np.sin(rad)**2;root=np.sqrt(magic)
 dlat=dlat*180/((6378137*(1-.00669342162296594323))/(magic*root)*p)
 dlng=dlng*180/(6378137/root*np.cos(rad)*p)
 inside=(lon>=72.004)&(lon<=137.8347)&(lat>=.8293)&(lat<=55.8271)
 result=c.copy();result[:,0]=lon+np.where(inside,dlng,0);result[:,1]=lat+np.where(inside,dlat,0)
 return result
def main():
    # Verify against independent scalar implementation at representative points.
    import eviltransform
    for lon,lat in [(121.4737,31.2304),(121.4365,31.195),(121.603,31.252)]:
     actual=convert(np.array([[lon,lat]]))[0];y,x=eviltransform.wgs2gcj(lat,lon)
     assert np.allclose(actual,[x,y],atol=1e-11,rtol=0)
    items=[('city_street_trade_areas.geojson','trade_areas_gcj02.geojson'),('city_street_partition.geojson','city_partition_gcj02.geojson'),('city_seed_anchors.geojson','seed_anchors_gcj02.geojson'),('unresolved_large_blocks.geojson','unresolved_blocks_gcj02.geojson'),('city_passage_links.gpkg','passage_links_gcj02.geojson')]
    manifest={'coordinate_system':'GCJ-02 approximate', 'source_assumption':'Source coordinates are interpreted as WGS84; administrative boundary datum remains unverified.', 'axis_order':'longitude, latitude','note':'GCJ-02 is not EPSG:4326. These are non-standard GeoJSON coordinate values for compatible basemaps, not RFC 7946 WGS84.', 'layers':[]}
    for source,target in items:
     print('processing',source,flush=True);df=gpd.read_file(src/source).to_crs(4326);original_count=len(df);df.geometry=shapely.transform(df.geometry.array,convert,include_z=False)
     # Unset CRS: never assert GCJ-02 is EPSG:4326.
     df.crs=None
     with warnings.catch_warnings():
      warnings.filterwarnings('ignore',message='.*without a CRS.*');df.to_file(out/target,driver='GeoJSON')
     manifest['layers'].append({'source':source,'output':target,'feature_count':original_count,'geometry_types':df.geom_type.value_counts().to_dict()})
    (out/'manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2));print(json.dumps(manifest,ensure_ascii=False,indent=2))

if __name__ == "__main__":
    main()
