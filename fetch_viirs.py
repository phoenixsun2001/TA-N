from pathlib import Path
import rasterio
from rasterio.windows import from_bounds
url='https://zenodo.org/records/17057876/files/nightlights.average_viirs.v21_m_500m_s_20230101_20231231_go_epsg4326_v20250904.tif?download=1'
Path('imagery').mkdir(exist_ok=True)
with rasterio.Env(GDAL_DISABLE_READDIR_ON_OPEN='EMPTY_DIR',GDAL_HTTP_TIMEOUT='60'):
 with rasterio.open(url) as src:
  win=from_bounds(121.35,31.10,121.70,31.34,src.transform).round_offsets().round_lengths()
  a=src.read(1,window=win); profile=src.profile.copy(); profile.update(width=int(win.width),height=int(win.height),transform=src.window_transform(win),compress='deflate')
  with rasterio.open('imagery/viirs_2023_native.tif','w',**profile) as dst:
   dst.write(a,1); dst.update_tags(source=url,source_scale_note='10x transformed; use raw as relative index, not radiance')
