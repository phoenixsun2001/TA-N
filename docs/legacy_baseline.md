# 上海全域影像融合 PoC

已完成 VIIRS、Sentinel-1、Sentinel-2 全市候选分区与消融对照，详见 [全市影像 PoC 说明](outputs/city_imagery/README.md)。人工与业务准确性验证尚未完成。

# 上海影像融合 PoC 更新

已新增实际使用 VIIRS、Sentinel-1 和 Sentinel-2 的三个上海试点，详见 [影像融合 PoC 说明](outputs/imagery_poc/README.md)。以下为此前的全市 POI 基线，并非影像融合成果。

# 上海自然商圈 PoC：公开数据基线

本项目以公开 OSM 数据和上海区界跑通一个**候选分区基线**。可用于检验 POI 聚集、阻隔近似及全域拓扑；结果不是已验证的自然商圈，也不能直接用于选址决策。对应原方案第一阶段的简化实现，不含真实步行网络、跨障碍通道或遥感融合。

## 查看结果

- `outputs/overview.png`：全市概览。彩色为候选归属，红线为候选核心，浅灰为未归属，深灰为阻隔代理，蓝色为水域代理。
- `outputs/commercial_core.geojson`：商业核心候选面；`id` 与中心一致。
- `outputs/trade_area_partition.geojson`：覆盖输入边界的分类面；正数 ID 为候选中心的归属区域，负数分别为未归属（-1）、阻隔代理（-2）、水域代理（-3）。
- `outputs/centers.geojson`：候选中心；`outputs/commercial_poi.geojson`：实际参与计算的 POI。
- `outputs/metrics.json`：本次数据校验、参数和面积。GeoJSON 输出为 WGS84；计算使用 UTM 51N 米制投影。

使用 QGIS 打开图层，按 `zone_type` 或 `id` 分类着色。中心 ID 仅是本次运行的内部标识，不是稳定的商圈名称或长期 ID。

## 算法及参数

1. 输入边界采用阿里云 DataV 地区边界的 16 个区的并集；输入 OSM 点中筛选餐饮、购物、休闲及部分生活服务类别。没有对商业 AOI 与点做合并去重，也没有闭店验证。
2. 将 POI 计数放到 250 米格网，做 sigma=1.2 格的高斯平滑（约 300 米标准差）；平滑结果 ≥0.40 且局部极大值为候选中心，中心间至少 1 公里。该阈值和间距只是 PoC 参数，没有按城乡分层校准。
3. 对地面高速、快速路和铁路进行栅格化作为阻隔代理；水域面单独栅格化。道路、铁路采取 all-touched 规则，所以可能放大阻隔。桥、隧道及非零层级本身不被视作地面阻隔；**没有建立桥梁、过街或地下通道连接**。这会使部分实际连通区域错误分离，网格精度也会使狭窄障碍丢失或变宽。
4. 从中心使用四邻域多源等成本扩张，最大累计格网距离 2 公里，其他格子为未归属；阻隔格子和水域格子另列为代理类型。该距离不是步行路网距离，不代表通行时间或消费者实际活动范围。核心区为归属内密度≥0.40的格子；不强求全域覆盖。
5. 将格网面裁切到边界并按 ID 合并，检查覆盖缺口、面重叠和几何有效性。计算结果存在多部件，因此同一 ID 的几何连续性也需人工核查。

这里的水域和阻隔是粗格网代理，并非精确水岸或道路多边形；输入区界是公开简化边界，不代表权威测绘边界。需要核实该边界的坐标基准、精度、授权与 OSM 坐标匹配，尤其注意中国地图常见坐标偏移。地块、机场、港区等功能区尚未识别，未归属不能被解释为无商业或无人口。OSM 记录并非完整商业清单，商业稀疏可能仅意味着数据缺失；结果不具备召回率或准确率评估。

## 数据来源

- [Geofabrik 上海 OSM 提取](https://download.geofabrik.de/asia/china/shanghai.html)，本地输入 `data/raw/shanghai.zip`；OpenStreetMap 数据 © OpenStreetMap contributors，使用须遵守 [ODbL](https://www.openstreetmap.org/copyright)。该提取的 POI 分类是 Geofabrik 简化图层，并非完整 OSM 标签。
- [DataV 上海区界](https://geo.datav.aliyun.com/areas_v3/bound/310000_full.json)，本地 `data/raw/shanghai_boundary.geojson`；用途和再发布许可须向数据提供方核实。
- 两份文件的 SHA-256 记录在 `outputs/metrics.json`，确保能核对具体输入快照。

## 重新运行

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
.venv/bin/python poc.py
```

已有输入文件放在 `data/raw/`，OSM 压缩包需解压至 `data/raw/osm/`。重新运行会覆盖 `outputs/` 文件。若需要重取输入：

```bash
mkdir -p data/raw/osm outputs
curl -L -o data/raw/shanghai.zip https://download.geofabrik.de/asia/china/shanghai-latest-free.shp.zip
unzip -o data/raw/shanghai.zip -d data/raw/osm
curl -L -o data/raw/shanghai_boundary.geojson https://geo.datav.aliyun.com/areas_v3/bound/310000_full.json
```

## 与目标方案的差距及验收

- 还没有 NPP-VIIRS、Sentinel-1、Sentinel-2：不能声称完成遥感融合，也没有夜光或地表指数。
- 尚无商业 AOI、可靠步行网络和出入口、障碍穿越通道，故无法评估真实步行可达性，边界未做道路级修整。
- 没有用徐家汇、陆家嘴、金桥三个试点的标注样本、订单或门店数据进行商业正确性验证。拓扑无缝不等于商业归属正确。
- 下一步应接入覆盖完整且授权可用的商业 POI / AOI、带跨越点的步行路网和核实坐标的行政边界；在三个试点人工标注拆分/合并与通道，再做阻隔违规抽查及覆盖率评估；再接入三类遥感质量控制和单元统计，比较融合前后误差。
