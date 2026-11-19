from pathlib import Path
import geopandas as g
from shapely.strtree import STRtree
out=Path('outputs/block_fit');blocks=g.read_file(out/'city_blocks.gpkg').to_crs(32651);anchors=g.read_file(out/'city_seed_anchors.gpkg').to_crs(32651);ix=STRtree(blocks.geometry.to_numpy());hits={}
for row in anchors[~anchors.matched].itertuples():
 for k in ix.query(row.geometry,predicate='intersects'):
  if blocks.area.iloc[k]>2e6:
   hits.setdefault(int(k),[]).append(int(row.seed_id));break
records=[{'block_id':int(blocks.id.iloc[k]),'area_km2':float(blocks.area.iloc[k]/1e6),'unmatched_count':len(v),'seed_ids':','.join(map(str,v)),'geometry':blocks.geometry.iloc[k]} for k,v in hits.items()]
g.GeoDataFrame(records,crs=32651).to_crs(4326).to_file(out/'unresolved_large_blocks.geojson',driver='GeoJSON');print('unresolved blocks',len(records),'seeds',sum(r['unmatched_count'] for r in records))
