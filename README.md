# 上海街区划分 PoC

使用 VIIRS、Sentinel-1、Sentinel-2，结合公开 POI、道路、铁路和水系，生成上海街区候选。当前交付为 v18 断头路延长切分实验版。

- [最终结果数据包（GCJ-02、原投影及审查图层）](outputs/deadend_extend_v18/shanghai_deadend_extension_gcj02.zip)
- [制作方法与质量说明](outputs/method_summary/上海街区制作方法.md)
- [流程图 Mermaid 源文件](outputs/method_summary/流程图.mmd)
- [早期基线说明](docs/legacy_baseline.md)

数据包包含原投影 EPSG:32651 GeoPackage 和近似 GCJ-02 GeoJSON。GCJ-02 用于显示，不应当作 WGS84 常规重投影。当前街区不是已验证的生活板块或商圈，存在既有面重叠、辅助切分线，以及 GCJ-02 和原投影范围不完全一致的问题，详见制作方法。

## 代码与数据

根目录 Python 脚本记录各阶段实验流程。可以使用 Python 3.9 和 `requirements.txt` 创建运行环境。脚本使用相对路径和阶段中间结果，不能仅凭此仓库直接从头运行完整流水线。

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
```

为控制仓库大小，本仓库包含最终版数据包和方法文档，不包含原始 OSM、用户提供的道路铁路水系、VIIRS/Sentinel 影像栅格、其他中间数据或 Python 虚拟环境。复现时需要按各脚本路径准备相应数据，并确认数据的许可和使用条件。公开 OSM 数据请遵守 [OpenStreetMap 版权与许可](https://www.openstreetmap.org/copyright)，行政边界来源与坐标精度仍需核实。
