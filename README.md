# uspex_parse
用于对uspex的结果进行解析后处理\
项目将使用8个碳的搜索结果进行处理
## 使用版本
USPEX_2021.0.1 (python版)\
Python 3.6.7 :: Anaconda custom (64-bit)
## results文件说明
BESTIndividuals --- 每一代最好的结构简略数据\
Individuals --- 所有生成结构的简略数据\
BESTgatheredPOSCARS --- 每一代最好的结构晶格坐标\
gatheredPOSCARS --- 所有的结构晶格坐标\
gatheredPOSCARS_unrelaxed ---所有的结构数据坐标（没做结构松弛）\
goodStructures --- 比较好的结构简略数据\
goodStructures_POSCARS ---比较好的结构晶格坐标\
origin ---每一代结构的来源\
parameters.uspex --- input.uspex源文件\
OUTPUT.txt --- 运行日志\
enthalpies_complete.csv --- 所有结构每个结构松弛阶段中的焓
### 若干矢量图
---
VarOperators.svg --- 每一代结构来源的占比\
E_series.svg --- 松弛步骤i和i+1的能量之间的相关性；有助于发现问题并改善结构松弛
enthalpy(raw)_vs_ID(raw).svg --- 焓-ID图\
enthalpy(per_atom)_vs_ID(raw).svg --- 每个原子的焓-ID图\
enthalpy(per_atom)_vs_cellUtility.volume(per_atom).svg --- 每个原子的焓-每个原子所占体积图\
enthalpy(per_atom)_statistics.svg --- 每个原子的焓统计图
## 参考文献
[[1]First principles prediction of amorphous phases using evolutionary algorithms](https://aip.scitation.org/doi/10.1063/1.4955105) \
[[2]A first principles study of amorphous and crystalline silicon tetraboride](https://www.sciencedirect.com/science/article/pii/S0254058420312888)