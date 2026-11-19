import json,numpy as np,geopandas as g,shapely
from shapely.strtree import STRtree
p='outputs/barrier_cut_v3/land_blocks_gcj02.geojson';d=g.read_file(p);a=d.geometry.to_numpy();print('read',len(a),flush=True);tree=STRtree(a);pairs=tree.query(a,predicate='overlaps').T;pairs=pairs[pairs[:,0]<pairs[:,1]];print('overlap pairs',len(pairs),flush=True)
areas=shapely.area(shapely.intersection(a[pairs[:,0]],a[pairs[:,1]])) if len(pairs) else np.array([])
result={'features':len(d),'invalid':int((~d.is_valid).sum()),'non_polygons':int((d.geom_type!='Polygon').sum()),'holes':int(sum(len(z.interiors)>0 for z in d.geometry)),'overlap_pairs':len(pairs),'overlap_area_degree2':float(areas.sum())};open('outputs/barrier_cut_v3/gcj_topology.json','w').write(json.dumps(result,indent=2));print(result,flush=True)
