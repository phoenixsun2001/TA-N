import io,urllib.request
import pyarrow.parquet as pq,pyarrow.compute as pc
import geopandas as g
from overturemaps.core import _create_s3_record_batch_reader
release='2026-10-21.0';bbox=(121.42,31.18,121.46,31.22);url=f'https://stac.overturemaps.org/{release}/collections.parquet';t=pq.read_table(io.BytesIO(urllib.request.urlopen(url,timeout=40).read()));x0,y0,x1,y1=bbox
flt=(pc.field('bbox','xmin')<x1)&(pc.field('bbox','xmax')>x0)&(pc.field('bbox','ymin')<y1)&(pc.field('bbox','ymax')>y0)
t=t.filter((pc.field('collection')=='segment')&(pc.field('type')=='Feature')&flt);paths=[a['aws']['alternate']['s3']['href'][5:] for a in t.column('assets').to_pylist()];print('indexed shards',len(paths),flush=True)
reader=_create_s3_record_batch_reader(paths,filter_expr=flt,connect_timeout=20,request_timeout=60)
if reader is None:raise RuntimeError('S3 access failed')
df=g.GeoDataFrame.from_arrow(reader);print('features',len(df),'columns',df.columns.tolist(),flush=True);df.to_parquet('data/raw/overture_xujiahui.parquet');print('saved')
