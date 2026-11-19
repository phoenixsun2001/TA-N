import sys
import geopandas as g,numpy as np,shapely
from shapely.strtree import STRtree
p=g.read_file(sys.argv[1] if len(sys.argv)>1 else 'outputs/revised_poc/activity_units_parent.gpkg').to_crs(32651)
# Count every exterior edge once after noding, and sample midpoints in <=25m pieces.
print('noding boundaries',flush=True);network=shapely.union_all(shapely.boundary(p.geometry.array));parts=shapely.get_parts(network);parts=shapely.segmentize(parts,25);coords,indices=shapely.get_coordinates(parts,return_index=True);valid=indices[1:]==indices[:-1];a=coords[:-1][valid];b=coords[1:][valid];mid=(a+b)*.5;length=np.linalg.norm(a-b,axis=1);valid=length>.001;mid=mid[valid];length=length[valid];points=shapely.points(mid);print('sample segments',len(length),'total km',length.sum()/1000,flush=True)
r=g.read_file('data/raw/osm/gis_osm_roads_free_1.shp').to_crs(32651);r=r[r.fclass.isin(['motorway','trunk','primary','secondary','tertiary','residential','unclassified','service','living_street'])&(r.bridge!='T')&(r.tunnel!='T')&(r.layer.fillna(0)==0)];
if len(sys.argv)>2:
 import pandas as pd
 extra=g.read_file(sys.argv[2]).to_crs(32651);r=g.GeoDataFrame(pd.concat([r,extra],ignore_index=True),crs=32651)
rail=g.read_file('data/raw/osm/gis_osm_railways_free_1.shp').to_crs(32651);rail=rail[rail.fclass.isin(['rail','light_rail','tram'])&(rail.bridge!='T')&(rail.tunnel!='T')&(rail.layer.fillna(0)==0)];water=g.read_file('data/raw/osm/gis_osm_water_a_free_1.shp').to_crs(32651);water=water[water.fclass.isin(['water','riverbank','reservoir','dock'])];city=g.read_file('data/raw/shanghai_boundary.geojson').to_crs(32651).geometry.union_all()
features=np.concatenate([r.geometry.to_numpy(),rail.geometry.to_numpy(),shapely.boundary(water.geometry.array),np.array([city.boundary],dtype=object)]);tree=STRtree(features);pairs=tree.query(points,predicate='dwithin',distance=3);covered=np.zeros(len(points),bool);covered[np.unique(pairs[0])]=True
# Strictly aligned with arbitrary 1 km raster cuts, excluding boundaries that are also near mapped physical features.
x=mid[:,0];y=mid[:,1];aligned=(((np.abs(x/1000-np.round(x/1000))*1000)<.05)&(np.abs(a[valid,0]-b[valid,0])<.05)) | (((np.abs(y/1000-np.round(y/1000))*1000)<.05)&(np.abs(a[valid,1]-b[valid,1])<.05))
roadtree=STRtree(r.geometry.to_numpy());hits=roadtree.query(points,predicate='dwithin',distance=3);roadmask=np.zeros(len(points),bool);roadmask[np.unique(hits[0])]=True
major=r[r.fclass.isin(['motorway','trunk','primary','secondary'])];mtree=STRtree(major.geometry.to_numpy());mhits=mtree.query(points,predicate='dwithin',distance=3);majormask=np.zeros(len(points),bool);majormask[np.unique(mhits[0])]=True
print('near physical (%)',round(length[covered].sum()/length.sum()*100,2),'near roads (%)',round(length[roadmask].sum()/length.sum()*100,2),'near major roads (%)',round(length[majormask].sum()/length.sum()*100,2),'not near km',round(length[~covered].sum()/1000,2),'grid aligned not near km',round(length[aligned&~covered].sum()/1000,2))
