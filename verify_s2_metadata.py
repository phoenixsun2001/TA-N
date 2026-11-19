import requests,re,json,planetary_computer as pc
from pystac_client import Client
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
path=Path('imagery/city_cache');m=json.loads((path/'sources.json').read_text());coll=Client.open('https://planetarycomputer.microsoft.com/api/stac/v1').get_collection('sentinel-2-l2a')
def f(x):
 item=coll.get_item(x['id']);r=requests.get(pc.sign(item.assets['product-metadata'].href),timeout=40);r.raise_for_status();t=r.text
 return {'id':x['id'],'offsets':sorted(set(re.findall(r'<BOA_ADD_OFFSET band_id="\d+">([^<]+)',t))),'quantification':re.findall(r'<BOA_QUANTIFICATION_VALUE[^>]*>([^<]+)',t)}
with ThreadPoolExecutor(max_workers=4) as p:data=list(p.map(f,m['sentinel2']))
(path/'s2_offsets.json').write_text(json.dumps(data,indent=2));print(json.dumps(data,indent=2))
