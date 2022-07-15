import os
import re
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from utils import *

def main():
    # 1.获取好的结构,并写入csv
    filename = './uspex/results1/goodStructures'
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
    db.to_csv('.\goodStructures.csv')

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

        if value <= bestEnthalpy + metastable:
            db3 = db3.append(db[db.index == index])

    db2.to_csv('.\stableStructures.csv')
    db3.to_csv('.\metastableStructures.csv')

def main1():
    # 1.先画一下焓分布图
    filename = '.\goodStructures.csv'
    db = pd.read_csv(filename)
    enthaply = np.array(db['Enthalpy (eV)'])
    enthaply /= 8
    print(enthaply)
    plt.bar(enthaply, np.ones_like(enthaply), fc = 'k', width = 1e-2)
    plt.xlabel('Energy/atom(eV)')
    plt.ylabel('Number Density')
    plt.yticks([])
    plt.savefig('./fig/enthalpy(per_atom)_statistics.png')
    plt.show()

if __name__ == '__main__':
    # main()
    # main1()
    # porcarSplit('')
    pos = readPOSCAR(r'.\POSCAR\EA18\POSCAR')
    print(pos)
    r, gr = radialDistributionFunction(pos)
    plt.plot(r, gr, lw = '1', c = 'r')
    plt.title('radial distribution function')
    plt.xlabel('r')
    plt.ylabel('g(r)')
    plt.show()
