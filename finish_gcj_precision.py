import geopandas as g,shapely,warnings,json
from shapely.geometry import Polygon
from shapely.strtree import STRtree
p='outputs/all_holes_gcj_v11/blocks_gcj02_final.geojson';d=g.read_file(p);a=d.geometry.to_numpy();tree=STRtree(a);out=[];filled=0;blocked=0
# Treat GCJ coordinates as a pseudo-WGS84 projection solely to measure numerical defects.
m=d.to_crs(32651);mg=m.geometry.to_numpy()
for i,z in enumerate(a):
 rings=[]
 for j,ring in enumerate(z.interiors):
  area=Polygon(mg[i].interiors[j]).area;hole=Polygon(ring);neighbors=tree.query(hole,predicate='intersects');occupied=sum(hole.intersection(a[k]).area for k in neighbors if k!=i)
  if area<=.01 and occupied<1e-17:filled+=1
  else:rings.append(ring);blocked+=1
 out.append(Polygon(z.exterior,rings))
d.geometry=out;d.crs=None
with warnings.catch_warnings():
 warnings.filterwarnings('ignore',message="'crs' was not provided.*");d.to_file('outputs/all_holes_gcj_v11/blocks_gcj02_clean.geojson',driver='GeoJSON')
print(json.dumps({'precision_holes_filled':filled,'holes_retained':blocked,'invalid':int((~d.is_valid).sum())},indent=2))
