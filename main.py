import os
import re
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from utils import *

def main():
    # 1.获取好的结构,并写入csv
    filename = r'.\uspex\results1\goodStructures'
    with open(filename, 'r') as f:
        text = f.read()
    text = text.replace(' ', '')
    # 正则表达式匹配出每一行数据
    pattern = re.compile(r"\|(\d+)\|\d+\|(\w+)\|(.*?)\|(.*?)\|(.*?)\|(.*?)\|(.*?)\|(.*?)\|(.*?)\|\n")
    match = pattern.findall(string = text)
    print(f'共有 {len(match)} 个较好结构')
    columns = ['ID', 'Origin', 'Composition', 'Enthalpy (eV)', 'Volume (A^3)', 'SYMMETRY (N)', 'Structure order', 'Average order', 'Quasientropy']
    db = pd.DataFrame(data = match, columns = columns)
    print(db.head())
    db.to_csv(r'.\goodStructures.csv', index = None)

    # 2.获取较稳定和亚稳定结构
    stable = 0.025852 # 1 / 2 k_b * t 在300k下为0.025852 eV
    metastable = 0.05
    db2 = pd.DataFrame()
    db3 = pd.DataFrame()
    bestEnthalpy = float(db['Enthalpy (eV)'][0])
    for index, value in enumerate(db['Enthalpy (eV)']):
        value = float(value)
        if value <= bestEnthalpy + stable:
            db2 = db2.append(db[db.index == index])

        if bestEnthalpy + stable <= value <= bestEnthalpy + metastable:
            db3 = db3.append(db[db.index == index])

    db2.to_csv(r'.\stableStructures.csv', index = None)
    db3.to_csv(r'.\metastableStructures.csv', index = None)

def main1():
    # 1.先画一下焓分布图
    filename = r'.\goodStructures.csv'
    db = pd.read_csv(filename)
    enthaply = np.array(db['Enthalpy (eV)'])
    enthaply /= 8
    print(enthaply)
    plt.bar(enthaply, np.ones_like(enthaply), fc = 'k', width = 1e-2)
    plt.xlabel('Energy/atom(eV)')
    plt.ylabel('Number Density')
    plt.yticks([])
    plt.savefig(r'.\fig\enthalpy(per_atom)_statistics.png')
    plt.show()
    # 2.划分一下POSCAR文件
    POSCARSplit(r'.\uspex\results1\goodStructures_POSCARS')

def main2():
    # 1.把每个结构的r和g(r)算出来先
    rs = []
    grs = []
    db = pd.read_csv('stableStructures.csv', index_col = None)
    dirs = os.listdir(r'.\POSCAR2')
    for dir in dirs:
        if int(dir[2:]) in list(db['ID']):
            pos = readPOSCAR(rf'.\POSCAR2\{dir}\POSCAR')
            r, gr = radialDistributionFunction(pos)
            rs.append(r)
            grs.append(gr)
    # 2.算出平均值
    r = max(rs, key = lambda x: len(x))
    gr = np.zeros_like(max(grs, key = lambda x: len(x)))
    count = np.zeros_like(gr)
    for g in grs:
        for index, value in enumerate(g):
            gr[index] += value
            count[index] += 1
    gr = np.divide(gr, count)
    drawRDF((r, gr, 'all_RDF'))


if __name__ == '__main__':
    # main()
    # main1()
    # main2()
    pass
    # drawRDFs(r'.\c8', False)
    # drawBADFs(r'.\c8', False)
    # print(distance(np.array([0.74775,0.81649,4.59397]), np.array([1.50980,0.57240,3.47602])))
    # print(np.arccos(np.linspace(-1,1)))
    drawBADF(arg = [r'.\dir_0.01\UCell_2_1.vasp'], isPrint = True)