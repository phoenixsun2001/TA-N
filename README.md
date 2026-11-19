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

## 最终版脚本

代码保留 v18 结果链及其影像、POI 和补充道路准备步骤。已有数据和中间结果准备齐全时，主要阶段按以下顺序执行：

1. `rebuild_road_level_v12.py`、`enrich_road_level_v12.py`：按道路等级生成初始面并关联 POI、影像。
2. `fill_all_holes_v13.py`、`close_deadend_notches_v14.py`：外环填充和内部狭长缺口处理。
3. `split_riverfront.py`：沿水细长区域候选切分。
4. `poi_anchored_cut.py`、`audit_poi_connections.py`、`split_direct_poi_v17.py`：码头和入口候选提取及直接切分。
5. `extend_deadends_v18.py`、`extend_notch_deadends_v18.py`：道路末端延长和边界凹口道路延长切分。

`convert_gcj02.py` 提供坐标转换；`poc.py`、`fetch_city.py`、`fetch_viirs.py`、`verify_s2_metadata.py`、`city_fusion.py` 和 `evaluate_supplement.py` 是此前影像与补充数据准备中使用的脚本。上列脚本依赖原始数据及部分历史生成的中间文件，例如影像栅格、筛选后的补充道路和人行通道结果；这些文件不在本仓库中，也没有完整的自动化重建入口。详见方法文档及各脚本输入路径。
