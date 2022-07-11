# uspex_parse

用于对uspex的结果进行解析后处理

项目将使用8个碳的搜索结果进行处理

## 使用版本
USPEX_2021.0.1 (python版)

## results文件说明
BESTIndividuals --- 每一代最好的结构简略数据

BESTgatheredPOSCARS --- 每一代最好的结构晶格坐标

gatheredPOSCARS --- 所有的结构晶格坐标

gatheredPOSCARS_unrelaxed ---所有的结构数据坐标（没做结构松弛）

goodStructures --- 比较好的结构简略数据

goodStructures_POSCARS ---比较好的结构晶格坐标

origin ---每一代结构的来源

parameters.uspex --- input.uspex源文件

VarOperators.svg --- 每一代结构来源的占比

OUTPUT.txt --- 运行日志

### 未知

Individuals

E_series.svg

enthalpy(raw)_vs_ID(raw).svg

enthalpies_complete.csv

enthalpy(per_atom)_vs_ID(raw).svg

enthalpy(per_atom)_vs_cellUtility.volume(per_atom).svg

enthalpy(per_atom)_statistics.svg
